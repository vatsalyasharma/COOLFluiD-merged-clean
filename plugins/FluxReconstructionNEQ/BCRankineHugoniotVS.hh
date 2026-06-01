#ifndef COOLFluiD_FluxReconstructionNEQ_BCRankineHugoniotVS_hh
#define COOLFluiD_FluxReconstructionNEQ_BCRankineHugoniotVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/VectorialFunction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Shock-fitting boundary condition for the downstream side of a fitted shock.
 *
 * Sets the ghost state to the UPSTREAM (freestream) state directly.
 * The Riemann solver at the shock face then resolves the jump between
 * the post-shock interior state and the freestream ghost state,
 * producing the correct shock flux.
 *
 * Unlike Dirichlet (which mirrors: ghost = 2*target - interior),
 * this BC sets ghost = target directly, because the interior and ghost
 * states are on opposite sides of a discontinuity.
 *
 * Usage in .CFcase:
 *   FluxReconstruction.Data.BcTypes = ... RankineHugoniotVS ...
 *   FluxReconstruction.Data.BcNames = ... ShockDown ...
 *   FluxReconstruction.Data.ShockDown.Vars = x y
 *   FluxReconstruction.Data.ShockDown.Def  = 0.0001952 0.004956 5590. 0. 1833.
 *   FluxReconstruction.ShockDown.applyTRS = ShockDown
 *
 * @author COOLFluiD Shock Fitting Extension
 */
class BCRankineHugoniotVS : public BCStateComputer {

public:

  static void defineConfigOptions(Config::OptionList& options);

  BCRankineHugoniotVS(const std::string& name);

  ~BCRankineHugoniotVS();

  static std::string getClassName() { return "BCRankineHugoniotVS"; }

  void setup();

  void unsetup();

  void configure(Config::ConfigArgs& args);

  void computeGhostStates(const std::vector< Framework::State* >& intStates,
                          std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& normals,
                          const std::vector< RealVector >& coords);

  void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                             std::vector< std::vector< RealVector* > >& ghostGrads,
                             const std::vector< RealVector >& normals,
                             const std::vector< RealVector >& coords);

protected:

  /// physical model var set
  Common::SafePtr<Framework::ConvectiveVarSet> m_varSet;

  /// transformer from input to update variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_inputToUpdateVar;

  /// user-defined functions for upstream (freestream) state
  std::vector<std::string> m_functions;

  /// variable names for function parser
  std::vector<std::string> m_vars;

  /// vectorial function parser
  Framework::VectorialFunction m_vFunction;

  /// input variable set name
  std::string m_inputVarStr;

  /// parsed input state
  Framework::State* m_inputState;

  /// space-time coordinate vector
  RealVector m_spaceTime;

  /// dimensional state after transformation
  Framework::State* m_dimState;

}; // class BCRankineHugoniotVS

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionNEQ_BCRankineHugoniotVS_hh
