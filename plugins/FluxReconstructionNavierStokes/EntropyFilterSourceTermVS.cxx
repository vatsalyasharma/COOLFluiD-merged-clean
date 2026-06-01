//Vatsalya: new class (FAILED APPROACH)-added S=-kappa*(I-F)*u as source term to make entropy filter dissipation visible to Newton Jacobian; ZERO EFFECT because PhysicalityCom runs before source terms so (I-F)*u~0 on already-filtered states; state-modifying limiters and residual source terms are incompatible when targeting same modal content

#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/EntropyFilterSourceTermVS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<EntropyFilterSourceTermVS, FluxReconstructionSolverData, FluxReconstructionNavierStokesModule>
    EntropyFilterSourceTermVSProvider("EntropyFilterSourceTermVS");

//////////////////////////////////////////////////////////////////////////////

void EntropyFilterSourceTermVS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Kappa","Filter relaxation rate (DynamicOption). Default 1.0.");
  options.addConfigOption< CFreal >("ZetaMax","Maximum filter strength. Default 36.0.");
  options.addConfigOption< CFreal >("Threshold","Modal energy ratio below which no filtering applied. Default 1e-4.");
  options.addConfigOption< CFreal >("Saturation","Modal energy ratio at which full filtering applied. Default 0.1.");
  options.addConfigOption< bool >("Active","Enable/disable at runtime (DynamicOption). Default true.");
  options.addConfigOption< CFuint >("ShowRate","Log active cell count every N iterations. Default 50.");
}

//////////////////////////////////////////////////////////////////////////////

EntropyFilterSourceTermVS::EntropyFilterSourceTermVS(const std::string& name) :
  StdSourceTerm(name),
  m_vdm(),
  m_vdmInv(),
  m_modeDegrees(),
  m_modalCoeffs(),
  m_dissipWeights(),
  m_kappa(1.0),
  m_zetaMax(36.0),
  m_threshold(1e-4),
  m_saturation(0.1),
  m_active(true),
  m_showRate(50),
  m_nbDims(0),
  m_nActiveCells(0)
{
  addConfigOptionsTo(this);

  setParameter("Kappa", &m_kappa);
  setParameter("ZetaMax", &m_zetaMax);
  setParameter("Threshold", &m_threshold);
  setParameter("Saturation", &m_saturation);
  setParameter("Active", &m_active);
  setParameter("ShowRate", &m_showRate);
}

//////////////////////////////////////////////////////////////////////////////

EntropyFilterSourceTermVS::~EntropyFilterSourceTermVS()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilterSourceTermVS::setup()
{
  CFAUTOTRACE;

  StdSourceTerm::setup();

  m_nbDims = PhysicalModelStack::getActive()->getDim();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // get Vandermonde matrix and inverse
  m_vdm    = *(frLocalData[0]->getVandermondeMatrix());
  m_vdmInv = *(frLocalData[0]->getVandermondeMatrixInv());

  const CFuint N = m_nbrSolPnts;
  const CFuint polyOrder = static_cast<CFuint>(frLocalData[0]->getPolyOrder());

  // compute mode degrees (same algorithm as EntropyFilteringVS)
  m_modeDegrees.resize(N);
  {
    const CFuint N_simplex_2D = (polyOrder + 1) * (polyOrder + 2) / 2;
    CFuint idx = 0;

    if (m_nbDims == 1)
    {
      for (CFuint d = 0; idx < N; ++d)
        m_modeDegrees[idx++] = d;
    }
    else if (N == N_simplex_2D && m_nbDims == 2)
    {
      for (CFuint d = 0; idx < N; ++d)
        for (CFuint m = 0; m < d + 1 && idx < N; ++m)
          m_modeDegrees[idx++] = d;
    }
    else
    {
      for (CFuint L = 0; idx < N; ++L)
        for (CFuint m = 0; m < 2 * L + 1 && idx < N; ++m)
          m_modeDegrees[idx++] = L;
    }
  }

  // allocate work arrays
  m_modalCoeffs.resize(N, m_nbrEqs);
  m_dissipWeights.resize(N);

  // log setup info
  CFLog(INFO, "EntropyFilterSourceTermVS::setup() - N=" << N << ", P=" << polyOrder
              << ", kappa=" << m_kappa << ", zetaMax=" << m_zetaMax
              << ", threshold=" << m_threshold << ", saturation=" << m_saturation
              << ", nbrEqs=" << m_nbrEqs << "\n");
  CFLog(INFO, "EntropyFilterSourceTermVS: mode degrees = [");
  for (CFuint k = 0; k < N; ++k)
  {
    CFLog(INFO, m_modeDegrees[k]);
    if (k < N - 1) CFLog(INFO, ", ");
  }
  CFLog(INFO, "]\n");
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilterSourceTermVS::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilterSourceTermVS::addSourceTerm(RealVector& resUpdates)
{
  if (!m_active)
  {
    resUpdates = 0.0;
    return;
  }

  const CFuint N = m_nbrSolPnts;

  // --- Step 1: Modal decomposition for all equations ---
  // û(k, iEq) = sum_j V^{-1}(k,j) * state_j[iEq]
  for (CFuint k = 0; k < N; ++k)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      CFreal sum = 0.0;
      for (CFuint j = 0; j < N; ++j)
        sum += m_vdmInv(k, j) * (*((*m_cellStates)[j]))[iEq];
      m_modalCoeffs(k, iEq) = sum;
    }
  }

  // --- Step 2: Compute modal energy ratio ---
  // E_mean = sum_iEq( û(0, iEq)^2 )
  // E_high = sum_iEq( sum_{k>0} û(k, iEq)^2 )
  // eta = E_high / (E_high + E_mean)
  CFreal E_mean = 0.0;
  CFreal E_high = 0.0;
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    E_mean += m_modalCoeffs(0, iEq) * m_modalCoeffs(0, iEq);
    for (CFuint k = 1; k < N; ++k)
      E_high += m_modalCoeffs(k, iEq) * m_modalCoeffs(k, iEq);
  }

  const CFreal E_total = E_mean + E_high;
  if (E_total < 1e-30)
  {
    resUpdates = 0.0;
    return;
  }

  const CFreal eta = E_high / E_total;

  // --- Step 3: Determine filter strength zeta from smoothness indicator ---
  if (eta < m_threshold)
  {
    // smooth cell: no filtering needed
    resUpdates = 0.0;
    return;
  }

  // linear ramp from threshold to saturation, clamped to [0, zetaMax]
  CFreal zetaCell;
  if (eta >= m_saturation)
  {
    zetaCell = m_zetaMax;
  }
  else
  {
    zetaCell = m_zetaMax * (eta - m_threshold) / (m_saturation - m_threshold);
  }

  ++m_nActiveCells;

  // --- Step 4: Compute dissipation weights d_k = 1 - exp(-zeta * p_k^2) ---
  m_dissipWeights[0] = 0.0;  // mean mode always preserved
  for (CFuint k = 1; k < N; ++k)
  {
    const CFreal pk = static_cast<CFreal>(m_modeDegrees[k]);
    m_dissipWeights[k] = 1.0 - std::exp(-zetaCell * pk * pk);
  }

  // --- Step 5: Compute source = -kappa * (I-F) * u ---
  // (I-F)*u[iSol][iEq] = sum_k V(iSol,k) * d_k * û(k, iEq)
  for (CFuint iSol = 0; iSol < N; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      CFreal val = 0.0;
      for (CFuint k = 0; k < N; ++k)
        val += m_vdm(iSol, k) * m_dissipWeights[k] * m_modalCoeffs(k, iEq);

      resUpdates[m_nbrEqs * iSol + iEq] = -m_kappa * val;
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilterSourceTermVS::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();

  // per-element-type hook: reset active cell counter at start of each type pass
  // and log the count from the previous iteration's accumulation
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  if (iter % m_showRate == 0 && m_iElemType == 0)
  {
    CFLog(NOTICE, "EntropyFilterSourceTermVS: " << m_nActiveCells
                  << " active cells (kappa=" << m_kappa
                  << ", zetaMax=" << m_zetaMax
                  << ", threshold=" << m_threshold
                  << ", saturation=" << m_saturation << ")\n");
  }
  m_nActiveCells = 0;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
