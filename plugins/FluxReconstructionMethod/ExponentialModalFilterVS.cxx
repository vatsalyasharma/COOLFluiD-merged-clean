//Vatsalya: new class-GFR-style exponential modal filter sigma(k)=exp(-alpha*(k/N)^s); precomputes filter matrix F=V*diag(sigma)*V^{-1} in setup(); placed in LimiterCom slot (postProcessSolutionImpl); all parameters are DynamicOption for runtime tuning

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"
#include "MathTools/MatrixInverter.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/ExponentialModalFilterVS.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ExponentialModalFilterVS, FluxReconstructionSolverData, FluxReconstructionModule>
    ExponentialModalFilterVSProvider("ExponentialModalFilterVS");

//////////////////////////////////////////////////////////////////////////////

void ExponentialModalFilterVS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("FilterStrength_NEW",
    "Exponential damping coefficient alpha in exp(-alpha * eta^s). Default 36.0 (GFR convention).");

  options.addConfigOption< CFuint >("FilterOrder_NEW",
    "Sharpness exponent s in exp(-alpha * eta^s). Default 0 means 2*P.");

  options.addConfigOption< CFuint >("ShowRate_NEW",
    "Print filter statistics every N iterations. Default 50.");

  options.addConfigOption< bool >("Active_NEW",
    "Enable/disable the filter at runtime (DynamicOption). Default true.");
}

//////////////////////////////////////////////////////////////////////////////

ExponentialModalFilterVS::ExponentialModalFilterVS(const std::string& name) :
  FluxReconstructionSolverCom(name),
  m_cellBuilder(CFNULL),
  m_nbrEqs(),
  m_nbrSolPnts(),
  m_nbFiltered(),
  m_totalNbFiltered()
{
  addConfigOptionsTo(this);

  m_filterStrength = 36.0;
  setParameter("FilterStrength_NEW", &m_filterStrength);

  m_filterOrder = 0;  // 0 means auto = 2*P
  setParameter("FilterOrder_NEW", &m_filterOrder);

  m_showRate = 50;
  setParameter("ShowRate_NEW", &m_showRate);

  m_active = true;
  setParameter("Active_NEW", &m_active);
}

//////////////////////////////////////////////////////////////////////////////

ExponentialModalFilterVS::~ExponentialModalFilterVS()
{
}

//////////////////////////////////////////////////////////////////////////////

void ExponentialModalFilterVS::configure(Config::ConfigArgs& args)
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ExponentialModalFilterVS::setup()
{
  CFAUTOTRACE;

  // get number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the local FR data
  vector<FluxReconstructionElementData*>& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);

  // get number of solution points
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  // get polynomial order
  const CFuint order = static_cast<CFuint>(frLocalData[0]->getPolyOrder());

  // auto filter order: 2*P if not specified
  if (m_filterOrder == 0)
  {
    m_filterOrder = 2 * order;
    if (m_filterOrder == 0) m_filterOrder = 2;  // safety for P0
  }

  // get Vandermonde matrix and its inverse
  RealMatrix vdm    = *(frLocalData[0]->getVandermondeMatrix());
  RealMatrix vdmInv = *(frLocalData[0]->getVandermondeMatrixInv());

  const CFuint N = m_nbrSolPnts;

  // build the diagonal filter weights sigma_k
  // sigma_k = exp(-alpha * (k/N)^s)  for k = 0..N-1
  // sigma_0 = 1 always (preserves cell average)
  RealVector sigma(N);
  sigma[0] = 1.0;

  for (CFuint k = 1; k < N; ++k)
  {
    const CFreal eta = static_cast<CFreal>(k) / static_cast<CFreal>(N);
    sigma[k] = exp(-m_filterStrength * pow(eta, static_cast<CFreal>(m_filterOrder)));
  }

  // build filter matrix: F = V * diag(sigma) * V^{-1}
  // First compute diag(sigma) * V^{-1}
  RealMatrix sigmaVdmInv(N, N);
  for (CFuint i = 0; i < N; ++i)
  {
    for (CFuint j = 0; j < N; ++j)
    {
      sigmaVdmInv(i, j) = sigma[i] * vdmInv(i, j);
    }
  }

  // Then compute F = V * (diag(sigma) * V^{-1})
  m_filterMatrix.resize(N, N);
  m_filterMatrix = 0.0;
  for (CFuint i = 0; i < N; ++i)
  {
    for (CFuint j = 0; j < N; ++j)
    {
      for (CFuint k = 0; k < N; ++k)
      {
        m_filterMatrix(i, j) += vdm(i, k) * sigmaVdmInv(k, j);
      }
    }
  }

  // allocate temp storage
  m_origStates.resize(N, m_nbrEqs);

  // log filter info
  CFLog(NOTICE, "ExponentialModalFilterVS: N=" << N << " P=" << order
        << " strength=" << m_filterStrength
        << " order=" << m_filterOrder << "\n");
  CFLog(NOTICE, "ExponentialModalFilterVS: sigma = [");
  for (CFuint k = 0; k < N; ++k)
  {
    CFLog(NOTICE, sigma[k]);
    if (k < N-1) CFLog(NOTICE, ", ");
  }
  CFLog(NOTICE, "]\n");
}

//////////////////////////////////////////////////////////////////////////////

void ExponentialModalFilterVS::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void ExponentialModalFilterVS::execute()
{
  CFTRACEBEGIN;

  // check if active (DynamicOption: can be toggled via interactive param file)
  if (!m_active)
  {
    CFTRACEEND;
    return;
  }

  // get the elementTypeData
  SafePtr<vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  const CFuint nbrElemTypes = elemType->size();
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  m_nbFiltered = 0;
  m_totalNbFiltered = 0;

  const CFuint N = m_nbrSolPnts;

  // loop over element types
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      geoData.idx = elemIdx;
      GeometricEntity* cell = m_cellBuilder->buildGE();

      vector<State*>* cellStates = cell->getStates();

      // copy original states
      for (CFuint iSol = 0; iSol < N; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_origStates(iSol, iEq) = (*((*cellStates)[iSol]))[iEq];
        }
      }

      // apply filter: U_new[i] = sum_j FilterMatrix(i,j) * U_old[j]
      for (CFuint iSol = 0; iSol < N; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          CFreal filtered = 0.0;
          for (CFuint jSol = 0; jSol < N; ++jSol)
          {
            filtered += m_filterMatrix(iSol, jSol) * m_origStates(jSol, iEq);
          }
          (*((*cellStates)[iSol]))[iEq] = filtered;
        }
      }

      ++m_nbFiltered;

      m_cellBuilder->releaseGE();
    }
  }

  // MPI reduction for statistics
  const std::string nsp = this->getMethodData().getNamespace();

#ifdef CF_HAVE_MPI
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  PE::GetPE().setBarrier(nsp);
  const CFuint count = 1;
  MPI_Allreduce(&m_nbFiltered, &m_totalNbFiltered, count, MPI_UNSIGNED, MPI_SUM, comm);
#else
  m_totalNbFiltered = m_nbFiltered;
#endif

  if (PE::GetPE().GetRank(nsp) == 0 && m_showRate > 0 && iter % m_showRate == 0)
  {
    CFLog(NOTICE, "ExponentialModalFilterVS: filtered " << m_totalNbFiltered
          << " cells (strength=" << m_filterStrength
          << ", order=" << m_filterOrder << ")\n");
  }

  PE::GetPE().setBarrier(nsp);

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
