#ifndef COOLFluiD_FluxReconstructionNEQ_BCSubInletSFInternalVS_hh
#define COOLFluiD_FluxReconstructionNEQ_BCSubInletSFInternalVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/VectorialFunction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Characteristic-based subsonic inflow BC for the downstream side of a
 * shock-fitted internal boundary. Full PLATO thermodynamics, no gamma
 * assumption.
 *
 * Prescribes (4 incoming characteristics):
 *   - h0 (total enthalpy, conserved across shock) from freestream
 *   - Y_s (species mass fractions, frozen across shock) from freestream
 *   - u_t = 0 (tangential velocity in face-normal frame)
 *
 * Extrapolates (1 outgoing characteristic):
 *   - M (Mach number, backward acoustic Riemann invariant)
 *   - s (entropy, from interior — not conserved across shock)
 *
 * Config:
 *   Def = rho_N rho_N2 u v T   (freestream state in Rhoivt)
 *   Vars = x y
 *
 * See doc/BCSubInletSFInternalVS_Design.md for full documentation.
 *
 * @author COOLFluiD Shock Fitting Extension
 */
class BCSubInletSFInternalVS : public BCStateComputer {

public:

  static void defineConfigOptions(Config::OptionList& options);

  BCSubInletSFInternalVS(const std::string& name);

  ~BCSubInletSFInternalVS();

  static std::string getClassName() { return "BCSubInletSFInternalVS"; }

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

  /// PLATO thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

  /// Update variable set (NEQ Rhoivt)
  Common::SafePtr<Framework::ConvectiveVarSet> m_varSet;

  /// Physical data arrays for interior and ghost
  RealVector m_intPhysData;
  RealVector m_ghostPhysData;

  /// Number of species
  CFuint m_nbSpecies;

  /// Freestream state [rho_N, rho_N2, u, v, T] (dimensional)
  std::vector<CFreal> m_fsState;

  /// Freestream total enthalpy h0 (dimensional, computed from PLATO)
  CFreal m_h0Target;

  /// Freestream species mass fractions
  RealVector m_ysTarget;

  /// Freestream mixture gas constant
  CFreal m_RmixTarget;

  /// Species gas constants [J/(kg K)]
  RealVector m_Ri;

  /// Working arrays for PLATO calls
  RealVector m_ys;       ///< species mass fractions (working)
  RealVector m_hsTot;    ///< species total enthalpies
  RealVector m_sSpecies; ///< species entropies
  RealVector m_tVec;     ///< temperature vector for PLATO

  /// User-defined functions for freestream state (spatial dependence)
  std::vector<std::string> m_functions;
  std::vector<std::string> m_vars;
  Framework::VectorialFunction m_vFunction;

  /// Space-time coordinate vector
  RealVector m_spaceTime;

  /// Parsed input state
  Framework::State* m_inputState;

  /// Variable transformer (input → update)
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_inputToUpdateVar;

  /// Input variable set name
  std::string m_inputVarStr;

  /// Reference values for adimensionalization
  RealVector m_refData;

  /// Max Newton iterations for temperature recovery
  static const CFuint MAX_NEWTON_ITERS = 25;

  /// Newton convergence tolerance
  static const CFreal NEWTON_TOL;

}; // class BCSubInletSFInternalVS

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionNEQ_BCSubInletSFInternalVS_hh
