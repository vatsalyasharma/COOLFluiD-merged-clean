#include "Framework/MethodStrategyProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/ParserException.hh"

#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNEQ/FluxReconstructionNEQ.hh"
#include "FluxReconstructionNEQ/BCSubInletSFInternalVS.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

const CFreal BCSubInletSFInternalVS::NEWTON_TOL = 1.0e-10;

Framework::MethodStrategyProvider<
    BCSubInletSFInternalVS,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNEQModule >
  BCSubInletSFInternalVSProvider("SubInletSFInternalVS");

//////////////////////////////////////////////////////////////////////////////

void BCSubInletSFInternalVS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def",
    "Freestream state in Rhoivt [rho_s1, rho_s2, ..., u, v, T] (dimensional).");
  options.addConfigOption< std::string >("InputVar","Input variable set (default: same as UpdateVar).");
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletSFInternalVS::BCSubInletSFInternalVS(const std::string& name) :
  BCStateComputer(name),
  m_library(CFNULL),
  m_varSet(CFNULL),
  m_intPhysData(),
  m_ghostPhysData(),
  m_nbSpecies(0),
  m_h0Target(0.0),
  m_ysTarget(),
  m_RmixTarget(0.0),
  m_Ri(),
  m_ys(),
  m_hsTot(),
  m_sSpecies(),
  m_tVec(),
  m_spaceTime(),
  m_inputState(CFNULL)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_functions = vector<std::string>();
  setParameter("Def", &m_functions);

  m_vars = vector<std::string>();
  setParameter("Vars", &m_vars);

  m_inputVarStr = "";
  setParameter("InputVar", &m_inputVarStr);
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletSFInternalVS::~BCSubInletSFInternalVS()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletSFInternalVS::computeGhostStates(
    const vector< State* >& intStates,
    vector< State* >& ghostStates,
    const std::vector< RealVector >& normals,
    const std::vector< RealVector >& coords)
{
  const CFuint nbrStates = intStates.size();
  cf_assert(nbrStates == ghostStates.size());
  cf_assert(nbrStates == normals.size());

  // Reference values from model (Fix #5: use model's own reference data)
  SafePtr<EulerTerm> eulerTerm = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();
  const RealVector& refData = eulerTerm->getReferencePhysicalData();
  const CFreal refRho = refData[EulerTerm::RHO];
  const CFreal refV   = refData[EulerTerm::V];
  const CFreal refT   = refData[EulerTerm::T];
  const CFreal refP   = refData[EulerTerm::P];
  // Note: refH = refV*refV (enthalpy ref = velocity² ref)

  // Freestream velocity in dimensional [m/s] for u_t computation (Fix #2)
  const CFreal u_fs_dim = m_fsState[m_nbSpecies];
  const CFreal v_fs_dim = m_fsState[m_nbSpecies + 1];

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    State& intSol   = *intStates[iState];
    State& ghostSol = *ghostStates[iState];
    const RealVector& normal = normals[iState];

    const CFreal nx = normal[XX];
    const CFreal ny = normal[YY];

    // ===== Step 1: Interior properties via physical data =====
    m_varSet->computePhysicalData(intSol, m_intPhysData);

    const CFreal u_int = m_intPhysData[EulerTerm::VX];  // non-dim
    const CFreal v_int = m_intPhysData[EulerTerm::VY];
    const CFreal V_int = m_intPhysData[EulerTerm::V];
    const CFreal a_int = m_intPhysData[EulerTerm::A];
    const CFreal T_int = m_intPhysData[EulerTerm::T];
    const CFreal p_int = m_intPhysData[EulerTerm::P];
    const CFreal rho_int = m_intPhysData[EulerTerm::RHO];
    const CFreal H_int = m_intPhysData[EulerTerm::H]; // non-dim total enthalpy

    // Interior species via physical data (Fix #6: use physData, not raw state)
    const CFuint firstSpecies = eulerTerm->getFirstScalarVar(0);
    for (CFuint i = 0; i < m_nbSpecies; ++i) {
      m_ys[i] = m_intPhysData[firstSpecies + i]; // mass fractions from physData
    }

    // Face-normal/tangential decomposition (non-dimensional)
    const CFreal u_n_int = u_int * nx + v_int * ny;
    const CFreal u_t_int = -u_int * ny + v_int * nx;

    // Mach number (non-dim V / non-dim a = dim V / dim a = correct)
    const CFreal Mach_int = (a_int > 1e-30) ? V_int / a_int : 0.0;

    // Total enthalpy (dimensional) = H_nondim * refV²
    const CFreal h0_int_dim = H_int * refV * refV;

    // Interior entropy from PLATO (dimensional)
    CFdouble T_int_dim = T_int * refT;
    CFdouble p_int_dim = p_int * refP;
    m_library->setSpeciesFractions(m_ys);
    m_library->getSpeciesEntropy(T_int_dim, p_int_dim, m_sSpecies);
    CFreal s_int = 0.0;
    for (CFuint i = 0; i < m_nbSpecies; ++i) s_int += m_ys[i] * m_sSpecies[i];

    // ===== Step 2: Mirror and extrapolate =====

    // Mirror h0 (total enthalpy conserved across shock)
    const CFreal h0_ghost_dim = 2.0 * m_h0Target - h0_int_dim;

    // Fix #2: tangential velocity = freestream tangential, NOT zero.
    // R-H jump: u_t is continuous across shock. Upstream u_t = freestream u_t.
    // Freestream velocity projected onto face tangent (dimensional):
    const CFreal u_t_fs_dim = -u_fs_dim * ny + v_fs_dim * nx;
    // Non-dimensionalize:
    const CFreal u_t_fs = u_t_fs_dim / refV;
    // Mirror: face-average u_t should equal freestream u_t
    const CFreal u_t_ghost = 2.0 * u_t_fs - u_t_int;

    // Extrapolate outgoing characteristic
    const CFreal Mach_ghost = Mach_int;
    const CFreal s_ghost = s_int;

    // Fix #4: species — set directly to freestream, no mirror.
    // The mirror fails when interior Y differs greatly from target (dissociation).
    // Species are incoming characteristics; prescribe directly.
    for (CFuint i = 0; i < m_nbSpecies; ++i) {
      m_ys[i] = m_ysTarget[i];
    }

    // Ghost R_mix
    CFreal R_mix_ghost = 0.0;
    for (CFuint i = 0; i < m_nbSpecies; ++i) R_mix_ghost += m_ys[i] * m_Ri[i];

    // Fix #3: Newton accounts for tangential KE.
    // Total velocity² = u_n² + u_t² = (M*a)² + u_t²
    // So: h(T) + 0.5*(M*a)² = h0 - 0.5*u_t²
    // Newton solves for T from the NORMAL component of KE only.
    const CFreal u_t_ghost_dim = u_t_ghost * refV;
    const CFreal h0_normal_dim = h0_ghost_dim - 0.5 * u_t_ghost_dim * u_t_ghost_dim;

    // ===== Step 3: Newton for T_ghost =====
    // Solve: h(T, Y_ghost) + 0.5*(M_ghost * a(T))² = h0_normal_dim

    // Vatsalya: This may reduce conv speed, use const gamma to test
    m_library->setSpeciesFractions(m_ys);

    CFdouble T_trial = T_int_dim;
    CFdouble p_trial = p_int_dim;
    CFdouble T_ghost_dim = T_trial;
    bool newtonConverged = false;

    for (CFuint iter = 0; iter < MAX_NEWTON_ITERS; ++iter)
    {
      // Enthalpy at trial T
      m_library->getSpeciesTotEnthalpies(T_trial, m_tVec, p_trial, m_hsTot);
      CFreal h_trial = 0.0;
      for (CFuint i = 0; i < m_nbSpecies; ++i) h_trial += m_ys[i] * m_hsTot[i];

      // Fix #7: recompute rho consistently for each T
      CFdouble rho_trial = p_trial / (R_mix_ghost * T_trial);
      CFdouble gamma_trial = 0.0, a_trial = 0.0;
      m_library->frozenGammaAndSoundSpeed(T_trial, p_trial, rho_trial,
                                          gamma_trial, a_trial, CFNULL);

      CFreal V_n_trial = Mach_ghost * a_trial;
      CFreal h0n_trial = h_trial + 0.5 * V_n_trial * V_n_trial;
      CFreal residual = h0n_trial - h0_normal_dim;

      if (std::abs(residual) < NEWTON_TOL * (std::abs(h0_normal_dim) + 1.0)) {
        T_ghost_dim = T_trial;
        newtonConverged = true;
        break;
      }

      // Fix #7: FD derivative with consistent rho at perturbed T
      const CFdouble dT_fd = 1.0;
      CFdouble T_plus = T_trial + dT_fd;
      m_library->getSpeciesTotEnthalpies(T_plus, m_tVec, p_trial, m_hsTot);
      CFreal h_plus = 0.0;
      for (CFuint i = 0; i < m_nbSpecies; ++i) h_plus += m_ys[i] * m_hsTot[i];

      CFdouble rho_plus = p_trial / (R_mix_ghost * T_plus); // Fix #7: recompute rho
      CFdouble gamma_plus = 0.0, a_plus = 0.0;
      m_library->frozenGammaAndSoundSpeed(T_plus, p_trial, rho_plus,
                                          gamma_plus, a_plus, CFNULL);
      CFreal V_n_plus = Mach_ghost * a_plus;
      CFreal h0n_plus = h_plus + 0.5 * V_n_plus * V_n_plus;
      CFreal dh0dT = (h0n_plus - h0n_trial) / dT_fd;

      if (std::abs(dh0dT) < 1e-30) {
        T_ghost_dim = T_trial;
        break;
      }

      T_trial -= residual / dh0dT;
      T_trial = std::max(T_trial, (CFdouble)100.0);
      T_ghost_dim = T_trial;
    }

    // Fix #7: warn on non-convergence
    if (!newtonConverged) {
      CFLog(WARN, "BCSubInletSFInternalVS: Newton T did not converge at iState="
            << iState << " T_ghost=" << T_ghost_dim << "\n");
    }

    // ===== Step 4: p_ghost from entropy =====
    CFdouble p_ref_d = 1.0e5;
    m_library->setSpeciesFractions(m_ys);
    m_library->getSpeciesEntropy(T_ghost_dim, p_ref_d, m_sSpecies);
    CFreal s_ref = 0.0;
    for (CFuint i = 0; i < m_nbSpecies; ++i) s_ref += m_ys[i] * m_sSpecies[i];

    CFreal p_ghost_dim = 1.0e5 * std::exp((s_ref - s_ghost) / R_mix_ghost);
    p_ghost_dim = std::max(p_ghost_dim, 1.0);

    // ===== Step 5: Set ghost state =====
    CFdouble rho_ghost_dim_d = p_ghost_dim / (R_mix_ghost * T_ghost_dim);
    CFdouble p_ghost_d = p_ghost_dim;
    CFdouble gamma_ghost = 0.0, a_ghost_dim = 0.0;
    m_library->frozenGammaAndSoundSpeed(T_ghost_dim, p_ghost_d, rho_ghost_dim_d,
                                        gamma_ghost, a_ghost_dim, CFNULL);

    const CFreal V_n_ghost_dim = Mach_ghost * a_ghost_dim;
    // Normal direction: same sign as interior u_n (inflow into downstream domain)
    const CFreal sign_un = (u_n_int >= 0.0) ? 1.0 : -1.0;
    const CFreal u_n_ghost_dim = sign_un * std::abs(V_n_ghost_dim);

    // Global velocity (dimensional)
    const CFreal u_ghost_dim = u_n_ghost_dim * nx - u_t_ghost_dim * ny;
    const CFreal v_ghost_dim = u_n_ghost_dim * ny + u_t_ghost_dim * nx;

    // Total enthalpy for physical data (dimensional)
    CFdouble T_gd = T_ghost_dim;
    m_library->getSpeciesTotEnthalpies(T_gd, m_tVec, p_ghost_d, m_hsTot);
    CFreal h_ghost = 0.0;
    for (CFuint i = 0; i < m_nbSpecies; ++i) h_ghost += m_ys[i] * m_hsTot[i];
    const CFreal V2_ghost_dim = u_ghost_dim * u_ghost_dim + v_ghost_dim * v_ghost_dim;
    const CFreal H_ghost_dim = h_ghost + 0.5 * V2_ghost_dim;

    // Fix #5: non-dimensionalize with model's refData, not ad-hoc refRho
    m_ghostPhysData[EulerTerm::RHO]   = (CFreal)rho_ghost_dim_d / refRho;
    m_ghostPhysData[EulerTerm::P]     = p_ghost_dim / refP;
    m_ghostPhysData[EulerTerm::T]     = (CFreal)T_ghost_dim / refT;
    m_ghostPhysData[EulerTerm::VX]    = u_ghost_dim / refV;
    m_ghostPhysData[EulerTerm::VY]    = v_ghost_dim / refV;
    m_ghostPhysData[EulerTerm::V]     = std::sqrt(V2_ghost_dim) / refV;
    m_ghostPhysData[EulerTerm::A]     = (CFreal)a_ghost_dim / refV;
    m_ghostPhysData[EulerTerm::GAMMA] = (CFreal)gamma_ghost;
    m_ghostPhysData[EulerTerm::H]     = H_ghost_dim / (refV * refV);
    m_ghostPhysData[EulerTerm::E]     = (h_ghost - p_ghost_dim / (CFreal)rho_ghost_dim_d)
                                        / (refV * refV);

    // Species mass fractions in physical data
    for (CFuint i = 0; i < m_nbSpecies; ++i) {
      m_ghostPhysData[firstSpecies + i] = m_ys[i];
    }

    // Convert to state via varSet (Fix #6: proper variable-set path)
    m_varSet->computeStateFromPhysicalData(m_ghostPhysData, ghostSol);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletSFInternalVS::computeGhostGradients(
    const std::vector< std::vector< RealVector* > >& intGrads,
    std::vector< std::vector< RealVector* > >& ghostGrads,
    const std::vector< RealVector >& normals,
    const std::vector< RealVector >& coords)
{
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads > 0);

  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    const RealVector& normal = normals[iState];
    const CFuint nbrGradVars = intGrads[iState].size();

    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      const RealVector& varGradI = *intGrads[iState][iGradVar];
      RealVector& varGradG = *ghostGrads[iState][iGradVar];
      const CFreal nVarGrad = MathFunctions::innerProd(varGradI, normal);
      varGradG = varGradI - 2.0 * nVarGrad * normal;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletSFInternalVS::setup()
{
  CFAUTOTRACE;
  BCStateComputer::setup();
  m_needsSpatCoord = true;

  // Variable set
  m_varSet = getMethodData().getUpdateVar();

  // Fix #6: assert Rhoivt (this BC only supports Rhoivt update variables)
  const std::string updateVarStr = getMethodData().getUpdateVarStr();
  if (updateVarStr.find("Rhoivt") == std::string::npos) {
    CFLog(WARN, "BCSubInletSFInternalVS: designed for Rhoivt update variables, "
          "got '" << updateVarStr << "'. Results may be incorrect.\n");
  }

  // Physical data arrays
  SafePtr<EulerTerm> eulerTerm = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();
  eulerTerm->resizePhysicalData(m_intPhysData);
  eulerTerm->resizePhysicalData(m_ghostPhysData);

  // PLATO library
  m_library = PhysicalModelStack::getActive()->getImplementor()->
    template getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(m_library.isNotNull());

  // Number of species
  m_nbSpecies = eulerTerm->getNbScalarVars(0);

  // Allocate
  m_ys.resize(m_nbSpecies);
  m_ysTarget.resize(m_nbSpecies);
  m_hsTot.resize(m_nbSpecies);
  m_sSpecies.resize(m_nbSpecies);
  m_Ri.resize(m_nbSpecies);
  m_tVec.resize(1); // CNEQ

  // Species gas constants
  m_library->setRiGas(m_Ri);

  // Parse freestream at (0,0,t=0)
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  m_spaceTime.resize(nbDims + 1);
  m_spaceTime = 0.0;

  const CFuint nbEq = PhysicalModelStack::getActive()->getNbEq();
  RealVector fsInput(nbEq);
  m_vFunction.evaluate(m_spaceTime, fsInput);

  // Store dimensional freestream state for u_t computation (Fix #2)
  m_fsState.resize(nbEq);
  for (CFuint i = 0; i < nbEq; ++i) m_fsState[i] = fsInput[i];

  // Freestream: [rho_s1, rho_s2, u, v, T] (dimensional)
  CFreal rho_fs = 0.0;
  for (CFuint i = 0; i < m_nbSpecies; ++i) rho_fs += fsInput[i];
  for (CFuint i = 0; i < m_nbSpecies; ++i) m_ysTarget[i] = fsInput[i] / rho_fs;

  const CFreal u_fs = fsInput[m_nbSpecies];
  const CFreal v_fs = fsInput[m_nbSpecies + 1];
  const CFreal T_fs = fsInput[m_nbSpecies + 2];
  const CFreal V_fs = std::sqrt(u_fs * u_fs + v_fs * v_fs);

  m_RmixTarget = 0.0;
  for (CFuint i = 0; i < m_nbSpecies; ++i) m_RmixTarget += m_ysTarget[i] * m_Ri[i];
  const CFreal p_fs = rho_fs * m_RmixTarget * T_fs;

  // h0 from PLATO (dimensional)
  m_library->setSpeciesFractions(m_ysTarget);
  CFdouble T_fs_d = T_fs;
  CFdouble p_fs_d = p_fs;
  m_library->getSpeciesTotEnthalpies(T_fs_d, m_tVec, p_fs_d, m_hsTot);
  CFreal h_fs = 0.0;
  for (CFuint i = 0; i < m_nbSpecies; ++i) h_fs += m_ysTarget[i] * m_hsTot[i];
  m_h0Target = h_fs + 0.5 * V_fs * V_fs;

  // Reference data for diagnostics
  const RealVector& refData = eulerTerm->getReferencePhysicalData();

  CFout << "BCSubInletSFInternalVS::setup()\n";
  CFout << "  Freestream: rho=" << rho_fs << " u=" << u_fs << " v=" << v_fs
        << " T=" << T_fs << " p=" << p_fs << "\n";
  CFout << "  Y_s = [";
  for (CFuint i = 0; i < m_nbSpecies; ++i) CFout << m_ysTarget[i] << (i<m_nbSpecies-1?", ":"");
  CFout << "]\n";
  CFout << "  R_mix=" << m_RmixTarget << " J/(kg K)\n";
  CFout << "  h0_target=" << m_h0Target << " J/kg\n";
  CFout << "  refRho=" << refData[EulerTerm::RHO] << " refV=" << refData[EulerTerm::V]
        << " refT=" << refData[EulerTerm::T] << " refP=" << refData[EulerTerm::P] << "\n";

  // Variable transformer
  const CFuint maxSt = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_inputToUpdateVar->setup(maxSt);
  m_inputState = new State();
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletSFInternalVS::unsetup()
{
  CFAUTOTRACE;
  deletePtr(m_inputState);
  BCStateComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletSFInternalVS::configure(Config::ConfigArgs& args)
{
  BCStateComputer::configure(args);

  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance(
      SubSystemStatusStack::getCurrentName()).getNamespace(namespc);
  SafePtr<PhysicalModel> physModel =
      PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string updateVarStr = getMethodData().getUpdateVarStr();
  if (m_inputVarStr.empty()) m_inputVarStr = updateVarStr;

  std::string provider = VarSetTransformer::getProviderName(
      physModel->getNameImplementor(), m_inputVarStr, updateVarStr);
  m_inputToUpdateVar =
      Environment::Factory<VarSetTransformer>::getInstance()
          .getProvider(provider)->create(physModel->getImplementor());

  m_vFunction.setFunctions(m_functions);
  if (m_vars.size() == 2 && m_vars[0] == "x" && m_vars[1] == "y") {
    m_vars.push_back("t");
  }
  m_vFunction.setVariables(m_vars);
  try { m_vFunction.parse(); }
  catch (Common::ParserException& e) { CFout << e.what() << "\n"; throw; }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
