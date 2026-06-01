//Vatsalya: new class-entropy-based positivity filter adapted from Dzanic & Witherden (JCP 468, 2022) for TCNEQ/CNEQ with PLATO mixture entropy; placed in PhysicalityCom slot (modifies states not residuals); zero Jacobian changes needed; includes Dzanic compliance fixes EF1-EF8 and MPI support

#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/EntropyFilteringVS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<EntropyFilteringVS, FluxReconstructionSolverData, FluxReconstructionNavierStokesModule>
    EntropyFilteringVSFRProvider("EntropyFilteringVS");

//////////////////////////////////////////////////////////////////////////////

void EntropyFilteringVS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("MinDensity","Minimum allowable value for species density.");
  options.addConfigOption< CFreal >("MinTemperature","Minimum allowable value for temperature (T and Tv).");
  options.addConfigOption< CFreal >("EntropyTol","Entropy tolerance for constraint relaxation.");
  options.addConfigOption< bool >("Active","DynamicOption: enable/disable the entropy filter at runtime.");
  options.addConfigOption< bool >("CheckFluxPnts","Whether to check flux point states in addition to solution points.");
  options.addConfigOption< CFuint >("BisectionIters","Number of bisection iterations for filter strength. Default 20.");
  // ShowRate is already defined in BasePhysicality — do NOT re-register
}

//////////////////////////////////////////////////////////////////////////////

EntropyFilteringVS::EntropyFilteringVS(const std::string& name) :
  BasePhysicality(name),
  socket_minEntropy("minEntropy"),
  m_neighborIDs(),
  m_speciesEntropy(),
  m_library(CFNULL),
  m_eulerVarSetMS(CFNULL),
  m_nbSpecies(0),
  m_TID(0),
  m_TvID(0),
  m_nbDims(0),
  m_nbTempVars(0),
  m_molMasses(),
  m_Rgas(8.314462618),
  m_nbrCells(0),
  m_cellMinEntropy(),
  m_stencilMinEntropy(0.0)
{
  addConfigOptionsTo(this);

  m_minDensity = 1e-16;
  setParameter( "MinDensity", &m_minDensity );

  m_minTemperature = 0.01;
  setParameter( "MinTemperature", &m_minTemperature );

  m_eTol = 1e-6;
  setParameter( "EntropyTol", &m_eTol );

  m_active = true;
  setParameter( "Active", &m_active );

  m_checkFlxPnts = true;
  setParameter( "CheckFluxPnts", &m_checkFlxPnts );

  m_nBisectionIters = 20;
  setParameter( "BisectionIters", &m_nBisectionIters );

  // m_showrate inherited from BasePhysicality (default=1), override default
  m_showrate = 50;
}

//////////////////////////////////////////////////////////////////////////////

EntropyFilteringVS::~EntropyFilteringVS()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilteringVS::configure ( Config::ConfigArgs& args )
{
  BasePhysicality::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
EntropyFilteringVS::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_posPrev);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
EntropyFilteringVS::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_outputPP);
  result.push_back(&socket_minEntropy);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilteringVS::setup()
{
  CFAUTOTRACE;

  // setup parent class
  BasePhysicality::setup();

  m_nbDims = PhysicalModelStack::getActive()->getDim();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // get multi-scalar term for species count
  m_eulerVarSetMS = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo< MultiScalarTerm< EulerTerm > >();
  cf_assert(!m_eulerVarSetMS.isNull());

  m_nbSpecies = m_eulerVarSetMS->getNbScalarVars(0);

  // detect variable set
  const std::string updateVarStr = getMethodData().getUpdateVarStr();
  if (updateVarStr == "RhoivtTv") {
    m_nbTempVars = 2;
  } else if (updateVarStr == "Rhoivt") {
    m_nbTempVars = 1;
  } else {
    CFLog(ERROR, "EntropyFilteringVS: unsupported update variable set '" << updateVarStr
                 << "'. Only RhoivtTv and Rhoivt are supported.\n");
    cf_assert(false);
  }

  m_TID  = m_nbSpecies + m_nbDims;
  m_TvID = m_nbSpecies + m_nbDims + 1;  // only used for RhoivtTv

  // get physical-chemical library
  m_library = PhysicalModelStack::getActive()->
    getImplementor()->template getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(!m_library.isNull());

  // get molar masses
  m_molMasses.resize(m_nbSpecies);
  m_library->getMolarMasses(m_molMasses);

  // resize working vectors
  m_speciesEntropy.resize(m_nbSpecies);

  // get the element type data
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();
  m_nbrCells = (*elemType)[0].getEndIdx();

  // initialize minEntropy socket (state-based for WriteSolutionHighOrder compatibility)
  DataHandle< CFreal > minEntropy = socket_minEntropy.getDataHandle();
  minEntropy.resize(m_nbrCells * m_nbrSolPnts);
  minEntropy = 0.0;

  // internal per-cell work array for neighbor exchange logic
  m_cellMinEntropy.resize(m_nbrCells, MathConsts::CFrealMax());

  // build neighbor IDs via face-adjacency (Dzanic 2022 Eq. 21)
  // Face-adjacent is the correct domain of influence for FR (not node-sharing)
  m_neighborIDs.resize(m_nbrCells);

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // use the face builder to iterate over interior faces
  SafePtr< GeometricEntityPool< FaceToCellGEBuilder > >
    faceBuilder = getMethodData().getFaceBuilder();
  FaceToCellGEBuilder::GeoData& geoDataFace = faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;

  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size() - 1;

  for (CFuint orient = 0; orient < nbrFaceOrients; ++orient)
  {
    const CFuint faceStartIdx = innerFacesStartIdxs[orient];
    const CFuint faceStopIdx  = innerFacesStartIdxs[orient + 1];

    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      geoDataFace.idx = faceID;
      GeometricEntity* face = faceBuilder->buildGE();

      const CFuint cellIDL = face->getNeighborGeo(LEFT )->getID();
      const CFuint cellIDR = face->getNeighborGeo(RIGHT)->getID();

      // Guard against ghost cells (ID >= m_nbrCells in MPI).
      // Only add local-to-local neighbors here; ghost cell entropy
      // is handled separately in execute() Pass 1b.
      const bool localL = (cellIDL < m_nbrCells);
      const bool localR = (cellIDR < m_nbrCells);

      if (localL && localR)
      {
        m_neighborIDs[cellIDL].push_back(cellIDR);
        m_neighborIDs[cellIDR].push_back(cellIDL);
      }

      faceBuilder->releaseGE();
    }
  }

  // ---- Modal filter infrastructure (Dzanic 2022) ----

  // get Vandermonde matrix and inverse
  m_vdm    = *(frLocalData[0]->getVandermondeMatrix());
  m_vdmInv = *(frLocalData[0]->getVandermondeMatrixInv());

  const CFuint N = m_nbrSolPnts;
  const CFuint polyOrder = static_cast<CFuint>(frLocalData[0]->getPolyOrder());

  // V[i,0] must be constant across all nodes (constant mode evaluates to same value everywhere)
  // For Legendre/quads: V[i,0] = 1.0. For Jacobi/triangles: V[i,0] ≈ 0.354
  m_vdm0 = m_vdm(0, 0);
  for (CFuint i = 1; i < N; ++i)
  {
    cf_assert(std::abs(m_vdm(i, 0) - m_vdm0) < 1e-10);
  }
  CFLog(INFO, "EntropyFilteringVS: V[i,0] = " << m_vdm0 << "\n");

  // compute mode degrees for exponential filter
  // For 2D simplex (triangle): modes ordered by total degree d,
  //   with (d+1) modes per degree. Degrees: 0, 1, 1, 2, 2, 2, ...
  // For 2D tensor-product (quad): modes ordered by level L=max(ksi,eta),
  //   with (2L+1) modes per level. Degrees: 0, 1, 1, 1, 2, 2, 2, 2, 2, ...
  // For 1D: modes ordered by degree k. Degrees: 0, 1, 2, ...
  m_modeDegrees.resize(N);
  {
    const CFuint N_simplex_2D = (polyOrder + 1) * (polyOrder + 2) / 2;
    CFuint idx = 0;

    if (m_nbDims == 1)
    {
      // 1D: mode k has degree k
      for (CFuint d = 0; idx < N; ++d)
        m_modeDegrees[idx++] = d;
    }
    else if (N == N_simplex_2D && m_nbDims == 2)
    {
      // 2D simplex: (d+1) modes of degree d
      for (CFuint d = 0; idx < N; ++d)
        for (CFuint m = 0; m < d + 1 && idx < N; ++m)
          m_modeDegrees[idx++] = d;
    }
    else
    {
      // 2D quad or unknown: (2L+1) modes of level L
      for (CFuint L = 0; idx < N; ++L)
        for (CFuint m = 0; m < 2 * L + 1 && idx < N; ++m)
          m_modeDegrees[idx++] = L;
    }
  }

  // allocate working arrays for bisection
  m_origStates.resize(N, m_nbrEqs);
  m_modalCoeffs.resize(N, m_nbrEqs);
  m_filteredStates.resize(N, m_nbrEqs);
  m_tempState.resize(m_nbrEqs);
  m_filterWeights.resize(N);

  // allocate filtered flux point states for bisection constraint checking
  m_filteredFlxPntStates.resize(m_maxNbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
    m_filteredFlxPntStates[iFlx].resize(m_nbrEqs);

  // log filter info
  CFLog(INFO, "EntropyFilteringVS::setup() - " << m_nbSpecies << " species, "
              << m_nbTempVars << " temp vars, " << m_nbrCells << " cells, "
              << m_nbrSolPnts << " sol pts/cell, P=" << polyOrder
              << ", bisection iters=" << m_nBisectionIters << "\n");
  CFLog(INFO, "EntropyFilteringVS: mode degrees = [");
  for (CFuint k = 0; k < N; ++k)
  {
    CFLog(INFO, m_modeDegrees[k]);
    if (k < N - 1) CFLog(INFO, ", ");
  }
  CFLog(INFO, "]\n");
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilteringVS::unsetup()
{
  CFAUTOTRACE;
  BasePhysicality::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilteringVS::execute()
{
  CFTRACEBEGIN;

  if (!m_active) {
    CFTRACEEND;
    return;
  }

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  const CFuint nbrElemTypes = elemType->size();
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  DataHandle< CFreal > minEntropy = socket_minEntropy.getDataHandle();
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();

  // if there is artificial viscosity, reset the positivity preservation values
  const bool hasArtVisc = getMethodData().hasArtificialViscosity();
  if (hasArtVisc)
  {
    DataHandle< CFreal > posPrev = socket_posPrev.getDataHandle();
    posPrev = MathConsts::CFrealMax();
  }

  // variable to store the number of limits done
  m_nbLimits = 0;
  m_totalNbLimits = 0;
  m_nbAvLimits = 0;
  m_totalNbAvLimits = 0;

  // ============================
  // Pass 1: Filter + compute new entropy
  // ============================
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();

    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // skip ghost cells (owned by another MPI rank)
      if (!((*m_cellStates)[0]->isParUpdatable()))
      {
        m_cellBuilder->releaseGE();
        continue;
      }

      // extrapolate states to flux points
      computeFlxPntStates(m_cellStatesFlxPnt);

      // read stencil minimum entropy from internal array (from previous iteration)
      m_stencilMinEntropy = m_cellMinEntropy[elemIdx];

      // check physicality
      if (!checkPhysicality())
      {
        // enforce physicality
        enforcePhysicality();

        ++m_nbLimits;
      }

      // Compute new per-cell minimum entropy across all solution points
      CFreal cellMin = MathConsts::CFrealMax();
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        const CFreal s = computeMixtureEntropy(*((*m_cellStates)[iSol]));
        cellMin = std::min(cellMin, s);
      }
      // Store as the new minimum for this cell (internal work array)
      m_cellMinEntropy[elemIdx] = cellMin;

      // release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  // ============================
  // Pass 1b: MPI ghost entropy exchange
  // ============================
  // For partition faces (one local cell, one ghost cell), compute the ghost
  // cell's minimum entropy from its already-synced states and fold it into
  // the local cell's m_cellMinEntropy. No explicit MPI send/recv needed —
  // ghost cell states are synced automatically by the framework.
  {
    SafePtr<TopologicalRegionSet> innerCells = MeshDataStack::getActive()->getTrs("InnerCells");
    SafePtr<TopologicalRegionSet> innerFaces = MeshDataStack::getActive()->getTrs("InnerFaces");

    SafePtr< GeometricEntityPool< FaceToCellGEBuilder > >
      faceBuilder = getMethodData().getFaceBuilder();
    FaceToCellGEBuilder::GeoData& geoDataFace = faceBuilder->getDataGE();
    geoDataFace.cellsTRS = innerCells;
    geoDataFace.facesTRS = innerFaces;
    geoDataFace.isBoundary = false;

    vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();
    const CFuint nbrFaceOrients = innerFacesStartIdxs.size() - 1;

    for (CFuint orient = 0; orient < nbrFaceOrients; ++orient)
    {
      const CFuint faceStartIdx = innerFacesStartIdxs[orient];
      const CFuint faceStopIdx  = innerFacesStartIdxs[orient + 1];

      for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
      {
        geoDataFace.idx = faceID;
        GeometricEntity* face = faceBuilder->buildGE();

        const CFuint cellIDL = face->getNeighborGeo(LEFT )->getID();
        const CFuint cellIDR = face->getNeighborGeo(RIGHT)->getID();

        const bool localL = (cellIDL < m_nbrCells);
        const bool localR = (cellIDR < m_nbrCells);

        // only process partition faces (one local, one ghost)
        if (localL != localR)
        {
          const CFuint ghostSide  = localL ? RIGHT : LEFT;
          const CFuint localCellID = localL ? cellIDL : cellIDR;

          // compute ghost cell minimum entropy from its states
          vector<State*>* ghostStates = face->getNeighborGeo(ghostSide)->getStates();
          CFreal ghostMinEntropy = MathConsts::CFrealMax();
          for (CFuint iSol = 0; iSol < ghostStates->size(); ++iSol)
          {
            const CFreal s = computeMixtureEntropy(*((*ghostStates)[iSol]));
            ghostMinEntropy = std::min(ghostMinEntropy, s);
          }

          // fold ghost entropy into local cell's stencil minimum
          m_cellMinEntropy[localCellID] = std::min(m_cellMinEntropy[localCellID], ghostMinEntropy);
        }

        faceBuilder->releaseGE();
      }
    }
  }

  // ============================
  // Pass 2: Neighbor entropy exchange
  // ============================
  // m_cellMinEntropy[cell] = min(m_cellMinEntropy[cell], m_cellMinEntropy[neighbor])
  // Copy first to avoid read-write race
  std::vector<CFreal> minEntCopy(m_cellMinEntropy);

  for (CFuint elemIdx = 0; elemIdx < m_nbrCells; ++elemIdx)
  {
    CFreal localMin = minEntCopy[elemIdx];
    for (CFuint iNeighbor = 0; iNeighbor < m_neighborIDs[elemIdx].size(); ++iNeighbor)
    {
      localMin = std::min(localMin, minEntCopy[m_neighborIDs[elemIdx][iNeighbor]]);
    }
    m_cellMinEntropy[elemIdx] = localMin;
  }

  // ============================
  // Write cell entropy to state-based socket for visualization
  // ============================
  // Each solution point in a cell gets the same cell-level value
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();

    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();
      m_cellStates = m_cell->getStates();

      // skip ghost cells (owned by another MPI rank)
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          minEntropy[((*m_cellStates)[iSol])->getLocalID()] = m_cellMinEntropy[elemIdx];
        }
      }

      m_cellBuilder->releaseGE();
    }
  }

  // ============================
  // MPI reduction + logging
  // ============================
  const std::string nsp = this->getMethodData().getNamespace();

#ifdef CF_HAVE_MPI
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  PE::GetPE().setBarrier(nsp);
  const CFuint count = 1;
  MPI_Allreduce(&m_nbLimits, &m_totalNbLimits, count, MPI_UNSIGNED, MPI_SUM, comm);
  MPI_Allreduce(&m_nbAvLimits, &m_totalNbAvLimits, count, MPI_UNSIGNED, MPI_SUM, comm);
#endif

  if (PE::GetPE().GetRank(nsp) == 0 && iter % m_showrate == 0)
  {
    CFLog(NOTICE, "EntropyFilteringVS: cells filtered = " << m_totalNbLimits
                  << ", avg limited = " << m_totalNbAvLimits << "\n");
  }

  PE::GetPE().setBarrier(nsp);

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

bool EntropyFilteringVS::checkPhysicality()
{
  bool physical = true;
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();

  // initialize output to 0
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    output[((*m_cellStates)[iSol])->getLocalID()] = 0.0;
  }

  // --- Check solution points ---
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const State& st = *((*m_cellStates)[iSol]);

    // check for NaN/Inf
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      if (!cfFinite(st[iEq])) { physical = false; break; }
    }
    if (!physical) break;

    // species density check
    for (CFuint i = 0; i < m_nbSpecies; ++i)
    {
      if (st[i] < m_minDensity) { physical = false; break; }
    }
    if (!physical) break;

    // temperature check (T and Tv if present)
    for (CFuint i = m_TID; i < m_nbrEqs; ++i)
    {
      if (st[i] < m_minTemperature) { physical = false; break; }
    }
    if (!physical) break;

    // entropy check (skip first iteration when stencilMinEntropy = +inf)
    if (m_stencilMinEntropy < 1e30)
    {
      const CFreal s = computeMixtureEntropy(st);
      if (s < m_stencilMinEntropy - m_eTol) { physical = false; break; }
    }
  }

  // --- Check flux points (optional) ---
  if (physical && m_checkFlxPnts)
  {
    for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
    {
      // NaN/Inf check
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        if (!cfFinite(m_cellStatesFlxPnt[iFlx][iEq])) { physical = false; break; }
      }
      if (!physical) break;

      // species density
      for (CFuint i = 0; i < m_nbSpecies; ++i)
      {
        if (m_cellStatesFlxPnt[iFlx][i] < m_minDensity) { physical = false; break; }
      }
      if (!physical) break;

      // temperature
      for (CFuint i = m_TID; i < m_nbrEqs; ++i)
      {
        if (m_cellStatesFlxPnt[iFlx][i] < m_minTemperature) { physical = false; break; }
      }
      if (!physical) break;

      // entropy at flux points
      if (m_stencilMinEntropy < 1e30)
      {
        const CFreal s = computeMixtureEntropy(m_cellStatesFlxPnt[iFlx]);
        if (s < m_stencilMinEntropy - m_eTol) { physical = false; break; }
      }
    }
  }

  return physical;
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilteringVS::enforcePhysicality()
{
  // =========================================================================
  // Dzanic & Witherden (JCP 468, 2022) Algorithm 2: Adaptive Entropy Filtering
  //
  // Bisect over filter strength ζ to find minimum exponential modal damping
  // that satisfies ALL constraints simultaneously:
  //   σ_k(ζ) = exp(-ζ * p_k²),  where p_k = degree of mode k
  //   ũ = V * diag(σ) * V^{-1} * u
  //
  // Constraints: ρ_i ≥ ρ_min,  T ≥ T_min,  Tv ≥ Tv_min,  s ≥ s_min - ε_σ
  // =========================================================================

  DataHandle< CFreal > output = socket_outputPP.getDataHandle();
  const CFuint N = m_nbrSolPnts;

  // --- Step 1: Save original nodal states ---
  for (CFuint iSol = 0; iSol < N; ++iSol)
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      m_origStates(iSol, iEq) = (*((*m_cellStates)[iSol]))[iEq];

  // --- Step 2: Compute modal coefficients û = V^{-1} * u ---
  for (CFuint k = 0; k < N; ++k)
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      CFreal sum = 0.0;
      for (CFuint j = 0; j < N; ++j)
        sum += m_vdmInv(k, j) * m_origStates(j, iEq);
      m_modalCoeffs(k, iEq) = sum;
    }

  // --- Step 3: Check cell mean (mode 0) physicality ---
  // At maximum filter strength, all states → V[i,0] * û[0].
  // If V[i,0]*û[0] is unphysical, no filter strength can produce a physical state.
  // For quads (Legendre): m_vdm0 = 1.0, so state = û[0] directly.
  // For triangles (Jacobi): m_vdm0 ≈ 0.354, actual state = 0.354 * û[0].
  bool meanFixed = false;

  // NaN/Inf sanitization
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    if (!cfFinite(m_modalCoeffs(0, iEq)))
    {
      m_modalCoeffs(0, iEq) = 0.0;
      meanFixed = true;
    }
  }
  // species density: check actual state value V[i,0]*û[0] against bound
  for (CFuint i = 0; i < m_nbSpecies; ++i)
  {
    if (m_vdm0 * m_modalCoeffs(0, i) < m_minDensity)
    {
      m_modalCoeffs(0, i) = 1.1 * m_minDensity / m_vdm0;
      meanFixed = true;
    }
  }
  // temperature: check actual state value V[i,0]*û[0] against bound
  for (CFuint i = m_TID; i < m_nbrEqs; ++i)
  {
    if (m_vdm0 * m_modalCoeffs(0, i) < m_minTemperature)
    {
      m_modalCoeffs(0, i) = 1.1 * m_minTemperature / m_vdm0;
      meanFixed = true;
    }
  }

  if (meanFixed)
  {
    m_nbAvLimits += 1;

    // Zero all modes except the fixed mean → all states become V[i,0] * û[0]
    for (CFuint iSol = 0; iSol < N; ++iSol)
    {
      output[((*m_cellStates)[iSol])->getLocalID()] = -10000.0;
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        (*((*m_cellStates)[iSol]))[iEq] = m_vdm0 * m_modalCoeffs(0, iEq);
    }
    return;
  }

  // --- Step 4: Bisection over ζ ∈ [0, ζ_max] ---
  // ζ = 0 → no filtering (unfiltered solution, which violates constraints)
  // ζ = ζ_max → maximum filtering (essentially cell mean, which passes)
  const CFreal zetaMax = 35.0;  // exp(-35) ≈ 6e-16 (FP64 safe)
  CFreal zeta_lo = 0.0;
  CFreal zeta_hi = zetaMax;

  for (CFuint bisecIter = 0; bisecIter < m_nBisectionIters; ++bisecIter)
  {
    const CFreal zeta = 0.5 * (zeta_lo + zeta_hi);

    // Compute filter weights: σ_k = exp(-ζ * p_k²)
    m_filterWeights[0] = 1.0;  // constant mode always preserved
    for (CFuint k = 1; k < N; ++k)
    {
      const CFreal pk = static_cast<CFreal>(m_modeDegrees[k]);
      m_filterWeights[k] = std::exp(-zeta * pk * pk);
    }

    // Apply filter: ũ_i = Σ_k V(i,k) * σ_k * û_k  for each equation
    for (CFuint iSol = 0; iSol < N; ++iSol)
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        CFreal val = 0.0;
        for (CFuint k = 0; k < N; ++k)
          val += m_vdm(iSol, k) * m_filterWeights[k] * m_modalCoeffs(k, iEq);
        m_filteredStates(iSol, iEq) = val;
      }

    // Check constraints on filtered states at solution points
    bool satisfied = true;

    for (CFuint iSol = 0; iSol < N && satisfied; ++iSol)
    {
      // species density: ρ_i ≥ ρ_min
      for (CFuint i = 0; i < m_nbSpecies; ++i)
      {
        if (m_filteredStates(iSol, i) < m_minDensity) { satisfied = false; break; }
      }
      if (!satisfied) break;

      // temperature: T, Tv ≥ T_min
      for (CFuint i = m_TID; i < m_nbrEqs; ++i)
      {
        if (m_filteredStates(iSol, i) < m_minTemperature) { satisfied = false; break; }
      }
      if (!satisfied) break;

      // entropy: s ≥ s_min_stencil - ε_σ
      if (m_stencilMinEntropy < 1e30)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
          m_tempState[iEq] = m_filteredStates(iSol, iEq);

        const CFreal s = computeMixtureEntropy(m_tempState);
        if (s < m_stencilMinEntropy - m_eTol) { satisfied = false; break; }
      }
    }

    // Check constraints at flux points (interpolated from filtered sol pts)
    // COOLFluiD uses open Gauss-Legendre nodes — flux points are NOT solution points
    if (satisfied && m_checkFlxPnts)
    {
      // interpolate filtered states to flux points
      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
      {
        m_filteredFlxPntStates[iFlx] = 0.0;
        for (CFuint iSol = 0; iSol < N; ++iSol)
        {
          const CFreal coef = (*m_solPolyValsAtFlxPnts)[iFlx][iSol];
          for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
            m_filteredFlxPntStates[iFlx][iEq] += coef * m_filteredStates(iSol, iEq);
        }
      }

      // check constraints at flux points
      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts && satisfied; ++iFlx)
      {
        for (CFuint i = 0; i < m_nbSpecies; ++i)
        {
          if (m_filteredFlxPntStates[iFlx][i] < m_minDensity) { satisfied = false; break; }
        }
        if (!satisfied) break;

        for (CFuint i = m_TID; i < m_nbrEqs; ++i)
        {
          if (m_filteredFlxPntStates[iFlx][i] < m_minTemperature) { satisfied = false; break; }
        }
        if (!satisfied) break;

        if (m_stencilMinEntropy < 1e30)
        {
          const CFreal s = computeMixtureEntropy(m_filteredFlxPntStates[iFlx]);
          if (s < m_stencilMinEntropy - m_eTol) { satisfied = false; break; }
        }
      }
    }

    // Update bisection brackets
    if (satisfied)
      zeta_hi = zeta;  // constraints OK → try less filtering
    else
      zeta_lo = zeta;  // constraints violated → need more filtering
  }

  // --- Step 5: Apply final filter with ζ = zeta_hi (guaranteed physical) ---
  m_filterWeights[0] = 1.0;
  for (CFuint k = 1; k < N; ++k)
  {
    const CFreal pk = static_cast<CFreal>(m_modeDegrees[k]);
    m_filterWeights[k] = std::exp(-zeta_hi * pk * pk);
  }

  for (CFuint iSol = 0; iSol < N; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      CFreal val = 0.0;
      for (CFuint k = 0; k < N; ++k)
        val += m_vdm(iSol, k) * m_filterWeights[k] * m_modalCoeffs(k, iEq);
      (*((*m_cellStates)[iSol]))[iEq] = val;
    }
    // Record filter strength for visualization
    output[((*m_cellStates)[iSol])->getLocalID()] = zeta_hi;
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal EntropyFilteringVS::computeMixtureEntropy(const RealVector& state)
{
  // For RhoivtTv/Rhoivt: state = [rho_0, rho_1, ..., rho_{NS-1}, u, v, T, (Tv)]
  // All quantities are non-dimensional. We dimensionalize for PLATO.
  //
  // Mixture entropy: s_mix = sum(Y_i * s_i) [J/(kg*K)]
  // where s_i comes from PLATO's species_entropy (partition function based).

  // get reference data
  RealVector& refData = m_eulerVarSetMS->getReferencePhysicalData();
  const CFreal refRho = refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO];
  const CFreal refT   = refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];

  // compute total density (non-dimensional)
  CFreal rhoTotal = 0.0;
  for (CFuint i = 0; i < m_nbSpecies; ++i)
  {
    rhoTotal += std::max(state[i], 0.0);
  }
  if (rhoTotal < 1e-30) return -1e30;  // degenerate state

  // dimensional temperature
  const CFreal T_dim = state[m_TID] * refT;
  if (T_dim < 1.0) return -1e30;  // degenerate state

  // compute mixture R = sum(Y_i / M_i) * R_universal
  CFreal sumYoverM = 0.0;
  for (CFuint i = 0; i < m_nbSpecies; ++i)
  {
    const CFreal Yi = std::max(state[i], 0.0) / rhoTotal;
    sumYoverM += Yi / m_molMasses[i];
  }
  const CFreal Rmix = m_Rgas * sumYoverM;

  // dimensional pressure and density
  const CFreal rho_dim = rhoTotal * refRho;
  const CFreal p_dim = rho_dim * Rmix * T_dim;

  if (p_dim < 1.0) return -1e30;  // degenerate state

  // get species entropy from PLATO [J/(kg*K)]
  CFdouble T_for_plato = T_dim;
  CFdouble p_for_plato = p_dim;
  m_library->getSpeciesEntropy(T_for_plato, p_for_plato, m_speciesEntropy);

  // validate: NaN/Inf from PLATO would poison entire neighbor stencil
  for (CFuint i = 0; i < m_nbSpecies; ++i)
  {
    if (!cfFinite(m_speciesEntropy[i])) return -1e30;
  }

  // mixture entropy: s_mix = sum(Y_i * s_i)
  CFreal sMix = 0.0;
  for (CFuint i = 0; i < m_nbSpecies; ++i)
  {
    const CFreal Yi = std::max(state[i], 0.0) / rhoTotal;
    sMix += Yi * m_speciesEntropy[i];
  }

  return sMix;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
