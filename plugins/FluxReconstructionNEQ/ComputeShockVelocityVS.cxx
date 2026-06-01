#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionNEQ/FluxReconstructionNEQ.hh"
#include "FluxReconstructionNEQ/ComputeShockVelocityVS.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolver.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeShockVelocityVS, DataProcessingData, FluxReconstructionNEQModule>
  computeShockVelocityVSProvider("ComputeShockVelocityVS");

//////////////////////////////////////////////////////////////////////////////

void ComputeShockVelocityVS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFreal> >("FreestreamDef",
    "Freestream state [rho_s1, rho_s2, ..., u, v, T] in Rhoivt (dimensional).");
}

//////////////////////////////////////////////////////////////////////////////

ComputeShockVelocityVS::ComputeShockVelocityVS(const std::string& name) :
  DataProcessingCom(name),
  socket_shockVelocity("shockVelocity"),
  socket_massFluxImbalance("massFluxImbalance"),
  socket_states("states"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  m_frData(CFNULL),
  m_faceBuilder(CFNULL),
  m_cellBuilder(CFNULL),
  m_rhoFs(0.0),
  m_uFs(0.0),
  m_vFs(0.0),
  m_nbrSolPnts(0),
  m_nbrFaceFlxPnts(0),
  m_nbSpecies(0),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_faceFlxPntConn(CFNULL),
  m_faceMappedCoordDir(CFNULL),
  m_flxLocalCoords(CFNULL),
  m_orient(0)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_fsDef = std::vector<CFreal>();
  setParameter("FreestreamDef", &m_fsDef);
}

//////////////////////////////////////////////////////////////////////////////

ComputeShockVelocityVS::~ComputeShockVelocityVS()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<SafePtr<BaseDataSocketSource>> ComputeShockVelocityVS::providesSockets()
{
  std::vector<SafePtr<BaseDataSocketSource>> result;
  result.push_back(&socket_shockVelocity);
  result.push_back(&socket_massFluxImbalance);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<SafePtr<BaseDataSocketSink>> ComputeShockVelocityVS::needsSockets()
{
  std::vector<SafePtr<BaseDataSocketSink>> result;
  result.push_back(&socket_states);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeShockVelocityVS::configure(Config::ConfigArgs& args)
{
  DataProcessingCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ComputeShockVelocityVS::setup()
{
  CFAUTOTRACE;
  DataProcessingCom::setup();

  // Get FR solver data
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<FluxReconstructionSolver> fr = spaceMethod.d_castTo<FluxReconstructionSolver>();
  cf_assert(fr.isNotNull());
  m_frData = fr->getData();

  // Get builders
  m_faceBuilder = m_frData->getFaceBuilder();
  m_cellBuilder = m_frData->getCellBuilder();

  // Get FR element data
  vector<FluxReconstructionElementData*>& frLocalData = m_frData->getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  // Get flux point data for face extrapolation (Fix #1: use flux-point states, not cell average)
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDir();
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  m_nbrFaceFlxPnts = frLocalData[0]->getNbrOfFaceFlxPnts();

  // Allocate flux point state storage
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_cellStatesFlxPnt.resize(m_nbrFaceFlxPnts);
  for (CFuint i = 0; i < m_nbrFaceFlxPnts; ++i)
  {
    m_cellStatesFlxPnt[i].resize(nbEqs);
  }
  m_unitNormalFlxPnts.resize(m_nbrFaceFlxPnts);
  for (CFuint i = 0; i < m_nbrFaceFlxPnts; ++i)
  {
    m_unitNormalFlxPnts[i].resize(PhysicalModelStack::getActive()->getDim());
  }

  // Parse freestream state
  cf_assert(m_fsDef.size() == nbEqs);
  m_nbSpecies = nbEqs - 3;
  m_rhoFs = 0.0;
  for (CFuint i = 0; i < m_nbSpecies; ++i) m_rhoFs += m_fsDef[i];
  m_uFs = m_fsDef[m_nbSpecies];
  m_vFs = m_fsDef[m_nbSpecies + 1];

  CFout << "ComputeShockVelocityVS::setup()\n";
  CFout << "  Freestream: rho=" << m_rhoFs << " u=" << m_uFs << " v=" << m_vFs << "\n";
  CFout << "  nbrSolPnts=" << m_nbrSolPnts << " nbrFaceFlxPnts=" << m_nbrFaceFlxPnts
        << " nbSpecies=" << m_nbSpecies << "\n";

  // Resize sockets
  DataHandle<CFreal> shockVel = socket_shockVelocity.getDataHandle();
  DataHandle<CFreal> mfluxImb = socket_massFluxImbalance.getDataHandle();
  const CFuint nbrStates = MeshDataStack::getActive()->getNbStates();
  shockVel.resize(nbrStates);
  shockVel = 0.0;
  mfluxImb.resize(nbrStates);
  mfluxImb = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeShockVelocityVS::unsetup()
{
  CFAUTOTRACE;
  m_cellStatesFlxPnt.clear();
  m_unitNormalFlxPnts.clear();
  DataProcessingCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeShockVelocityVS::executeOnTrs()
{
  CFAUTOTRACE;

  // Zero out sockets
  DataHandle<CFreal> shockVel = socket_shockVelocity.getDataHandle();
  DataHandle<CFreal> mfluxImb = socket_massFluxImbalance.getDataHandle();
  shockVel = 0.0;
  mfluxImb = 0.0;

  SafePtr<TopologicalRegionSet> cellTrs =
    MeshDataStack::getActive()->getTrs("InnerCells");
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // Get boundary face start indices
  std::map<std::string, std::vector<std::vector<CFuint>>>& bndFacesStartIdxsAll =
    m_frData->getBndFacesStartIdxs();
  std::vector<std::vector<CFuint>>& bndFacesStartIdxs =
    bndFacesStartIdxsAll[faceTrs->getName()];
  const CFuint nbOrients = bndFacesStartIdxs[0].size() - 1;
  const CFuint nbTRs = faceTrs->getNbTRs();

  // Face Jacobian vector sizes
  DataHandle<vector<CFreal>> faceJacobVecSizeFaceFlxPnts =
    socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();

  // Set up builders
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;

  CellToFaceGEBuilder::GeoData& geoDataCB = m_cellBuilder->getDataGE();
  geoDataCB.trs = cellTrs;

  CFuint nFaces = 0;
  CFreal maxVs = 0.0;

  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient + 1];

      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        geoData.idx = faceID;
        GeometricEntity* face = m_faceBuilder->buildGE();

        GeometricEntity* intCell = face->getNeighborGeo(0);
        std::vector<State*>* cellStates = intCell->getStates();

        if ((*cellStates)[0]->isParUpdatable())
        {
          // --- Fix #1: Extrapolate states to face flux points ---
          // Compute face Jacobian vectors for unit normals
          vector<RealVector> faceJacobVecs =
            face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

          // Extrapolate solution to flux points and compute normals
          for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
          {
            // Reset flux point state
            m_cellStatesFlxPnt[iFlx] = 0.0;

            // Flux point index in the element's local numbering
            const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];

            // Extrapolate: sum over solution points
            for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
            {
              const CFreal coeff = (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol];
              for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
              {
                m_cellStatesFlxPnt[iFlx][iEq] += coeff * (*(*cellStates)[iSol])[iEq];
              }
            }

            // Unit normal at this flux point
            const CFreal faceJacobSize = faceJacobVecSizeFaceFlxPnts[face->getID()][iFlx];
            m_unitNormalFlxPnts[iFlx] = faceJacobVecs[iFlx] / faceJacobSize;
          }

          // --- Compute Vs averaged over face flux points ---
          CFreal sumVs = 0.0;
          CFreal sumImbal = 0.0;

          for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
          {
            const CFreal nx = m_unitNormalFlxPnts[iFlx][XX];
            const CFreal ny = m_unitNormalFlxPnts[iFlx][YY];

            // Interior state at this flux point
            CFreal rhoInt = 0.0;
            for (CFuint i = 0; i < m_nbSpecies; ++i) rhoInt += m_cellStatesFlxPnt[iFlx][i];
            const CFreal uInt = m_cellStatesFlxPnt[iFlx][m_nbSpecies];
            const CFreal vInt = m_cellStatesFlxPnt[iFlx][m_nbSpecies + 1];

            // Normal velocities
            const CFreal unFs  = m_uFs * nx + m_vFs * ny;
            const CFreal unInt = uInt * nx + vInt * ny;

            // Mass fluxes
            const CFreal mfluxFs  = m_rhoFs * unFs;
            const CFreal mfluxInt = rhoInt * unInt;

            // Shock velocity
            const CFreal drho = m_rhoFs - rhoInt;
            CFreal Vs = 0.0;
            if (std::abs(drho) > 1e-30)
            {
              Vs = (mfluxFs - mfluxInt) / drho;
            }

            CFreal imbal = 0.0;
            if (std::abs(mfluxFs) > 1e-30)
            {
              imbal = 100.0 * (mfluxFs - mfluxInt) / std::abs(mfluxFs);
            }

            sumVs += Vs;
            sumImbal += imbal;
          }

          // Average over flux points on this face
          const CFreal avgVs = sumVs / (CFreal)m_nbrFaceFlxPnts;
          const CFreal avgImbal = sumImbal / (CFreal)m_nbrFaceFlxPnts;

          // Store in socket for all states of this cell
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            const CFuint stateID = (*cellStates)[iSol]->getLocalID();
            shockVel[stateID] = avgVs;
            mfluxImb[stateID] = avgImbal;
          }

          maxVs = std::max(maxVs, std::abs(avgVs));
          ++nFaces;
        }

        m_faceBuilder->releaseGE();
      }
    }
  }

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  CFout << "ComputeShockVelocityVS [" << faceTrs->getName() << "] iter " << iter
        << ": " << nFaces << " faces, max|Vs|=" << maxVs << " m/s\n";
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
