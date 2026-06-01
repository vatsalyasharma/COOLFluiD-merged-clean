#include "Framework/MethodStrategyProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "FluxReconstructionNEQ/FluxReconstructionNEQ.hh"
#include "FluxReconstructionNEQ/BCRankineHugoniotVS.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCRankineHugoniotVS,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNEQModule >
  BCRankineHugoniotVSProvider("RankineHugoniotVS");

//////////////////////////////////////////////////////////////////////////////

void BCRankineHugoniotVS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the upstream (freestream) state.");
  options.addConfigOption< std::string >("InputVar","Input variable set (default: same as UpdateVar).");
}

//////////////////////////////////////////////////////////////////////////////

BCRankineHugoniotVS::BCRankineHugoniotVS(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_inputToUpdateVar(),
  m_inputState(),
  m_spaceTime(),
  m_dimState()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_functions = vector<std::string>();
  setParameter("Def",&m_functions);

  m_vars = vector<std::string>();
  setParameter("Vars",&m_vars);

  m_inputVarStr = "";
  setParameter("InputVar",&m_inputVarStr);
}

//////////////////////////////////////////////////////////////////////////////

BCRankineHugoniotVS::~BCRankineHugoniotVS()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCRankineHugoniotVS::computeGhostStates(
    const vector< State* >& intStates,
    vector< State* >& ghostStates,
    const std::vector< RealVector >& normals,
    const std::vector< RealVector >& coords)
{
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  CFreal time = subSysStatus->getCurrentTimeDim();

  const CFuint nbrStates = intStates.size();
  cf_assert(nbrStates == ghostStates.size());
  cf_assert(nbrStates == normals.size());

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    State& ghostSol = *ghostStates[iState];

    // Evaluate upstream (freestream) state from user-defined function
    for (CFuint i = 0; i < coords[iState].size(); ++i)
    {
      m_spaceTime[i] = coords[iState][i];
    }
    m_spaceTime[coords[iState].size()] = time;
    m_vFunction.evaluate(m_spaceTime, *m_inputState);

    // Transform to update variables
    m_dimState = m_inputToUpdateVar->transform(m_inputState);

    // Adimensionalize and store
    m_varSet->setAdimensionalValues(*m_dimState, ghostSol);

    // NO MIRROR: ghost = upstream state directly.
    // Unlike Dirichlet which does ghost = 2*target - interior,
    // we set ghost = target because the interior and ghost are
    // on opposite sides of a shock discontinuity.
    // The Riemann solver resolves the jump.
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCRankineHugoniotVS::computeGhostGradients(
    const std::vector< std::vector< RealVector* > >& intGrads,
    std::vector< std::vector< RealVector* > >& ghostGrads,
    const std::vector< RealVector >& normals,
    const std::vector< RealVector >& coords)
{
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());
  cf_assert(nbrStateGrads > 0);

  // Zero normal gradient: flip the normal component of interior gradient.
  // This is appropriate for a shock boundary where the upstream flow
  // is uniform (freestream) and has zero gradient.
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    const RealVector& normal = normals[iState];
    const CFuint nbrGradVars = intGrads[iState].size();
    cf_assert(nbrGradVars == ghostGrads[iState].size());

    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      const RealVector& varGradI = *intGrads[iState][iGradVar];
      RealVector& varGradG = *ghostGrads[iState][iGradVar];

      const CFreal nVarGrad = MathTools::MathFunctions::innerProd(varGradI, normal);
      varGradG = varGradI - 2.0 * nVarGrad * normal;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCRankineHugoniotVS::setup()
{
  CFAUTOTRACE;

  BCStateComputer::setup();

  // need spatial coordinates to evaluate freestream function
  m_needsSpatCoord = true;

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_inputToUpdateVar->setup(maxNbStatesInCell);

  m_inputState = new State();

  m_varSet = getMethodData().getUpdateVar();

  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  m_spaceTime.resize(nbDims + 1);
}

//////////////////////////////////////////////////////////////////////////////

void BCRankineHugoniotVS::unsetup()
{
  CFAUTOTRACE;
  deletePtr(m_inputState);
  BCStateComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void BCRankineHugoniotVS::configure(Config::ConfigArgs& args)
{
  BCStateComputer::configure(args);

  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance(
      SubSystemStatusStack::getCurrentName()).getNamespace(namespc);
  SafePtr<PhysicalModel> physModel =
      PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string updateVarStr = getMethodData().getUpdateVarStr();

  if (m_inputVarStr.empty())
  {
    m_inputVarStr = updateVarStr;
  }

  std::string provider = VarSetTransformer::getProviderName(
      physModel->getNameImplementor(), m_inputVarStr, updateVarStr);

  m_inputToUpdateVar =
      Environment::Factory<VarSetTransformer>::getInstance()
          .getProvider(provider)->create(physModel->getImplementor());
  cf_assert(m_inputToUpdateVar.isNotNull());

  m_vFunction.setFunctions(m_functions);

  if (m_vars.size() == 2 && m_vars[0] == "x" && m_vars[1] == "y") {
    m_vars.push_back("t");
  }
  if (m_vars.size() == 3 && m_vars[0] == "x" && m_vars[1] == "y" && m_vars[2] == "z") {
    m_vars.push_back("t");
  }
  m_vFunction.setVariables(m_vars);

  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw;
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
