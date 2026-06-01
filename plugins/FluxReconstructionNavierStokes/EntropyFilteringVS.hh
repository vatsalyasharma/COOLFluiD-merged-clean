#ifndef COOLFluiD_FluxReconstructionMethod_EntropyFilteringVS_hh
#define COOLFluiD_FluxReconstructionMethod_EntropyFilteringVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/MultiScalarTerm.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/BasePhysicality.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class EulerTerm;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Entropy-based positivity filter for TCNEQ/CNEQ (RhoivtTv/Rhoivt).
 * Adapted from Dzanic & Witherden (JCP 468, 2022) for multi-species
 * thermochemical nonequilibrium using PLATO for mixture entropy.
 *
 * Enforces: species density > 0, T/Tv > T_min, mixture entropy >= stencil minimum.
 * Blends toward cell average when constraints are violated.
 *
 * @author Vatsalya Sharma
 */
class EntropyFilteringVS : public BasePhysicality {
public:

  /**
   * Constructor.
   */
  explicit EntropyFilteringVS(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~EntropyFilteringVS();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Setup private data
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Execute: two-pass algorithm (filter + neighbor exchange)
   */
  void execute();

  /// Returns the DataSocket's that this command provides as sources
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected: // functions

  /**
   * Check if the states are physical (density, temperature, entropy)
   */
  virtual bool checkPhysicality();

  /**
   * Enforce physicality via sequential constraint blending
   */
  virtual void enforcePhysicality();

  /**
   * Compute mixture entropy for a given state (RhoivtTv/Rhoivt)
   * @param state the solution state
   * @return mixture entropy s_mix = sum(Y_i * s_i) [J/(kg*K)]
   */
  CFreal computeMixtureEntropy(const RealVector& state);

protected: // data

  /// minimum allowable value for species density
  CFreal m_minDensity;

  /// minimum allowable value for temperature (T and Tv)
  CFreal m_minTemperature;

  /// entropy tolerance for constraint relaxation
  CFreal m_eTol;

  /// DynamicOption: active flag (can be toggled at runtime)
  bool m_active;

  /// whether to also check flux point states
  bool m_checkFlxPnts;

  /// number of bisection iterations for filter strength (default 20)
  CFuint m_nBisectionIters;

  /// socket for minimum entropy per cell (source)
  Framework::DataSocketSource< CFreal > socket_minEntropy;

  /// neighbor cell IDs for each cell (face-adjacent, Dzanic 2022 Eq. 21)
  std::vector< std::vector<CFuint> > m_neighborIDs;

  /// species entropies per unit mass [J/(kg*K)] (from PLATO)
  RealVector m_speciesEntropy;

  /// physical-chemical library handle (PLATO entropy)
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

  /// multi-scalar convective term
  Common::SafePtr< Framework::MultiScalarTerm< Physics::NavierStokes::EulerTerm > > m_eulerVarSetMS;

  /// number of species
  CFuint m_nbSpecies;

  /// index of T in state vector (= nbSpecies + dim)
  CFuint m_TID;

  /// index of Tv in state vector (= nbSpecies + dim + 1), only for RhoivtTv
  CFuint m_TvID;

  /// number of dimensions
  CFuint m_nbDims;

  /// number of temperature variables (1 for Rhoivt, 2 for RhoivtTv)
  CFuint m_nbTempVars;

  /// species molar masses [kg/mol]
  RealVector m_molMasses;

  /// universal gas constant [J/(mol*K)]
  CFreal m_Rgas;

  /// number of cells
  CFuint m_nbrCells;

  /// minimum entropy from previous iteration (per cell, internal work array)
  std::vector<CFreal> m_cellMinEntropy;

  /// stencil minimum entropy for current cell
  CFreal m_stencilMinEntropy;

  // ---- Modal filter infrastructure (Dzanic 2022, Algorithm 2) ----

  /// Vandermonde matrix V: modes → nodes [N × N]
  RealMatrix m_vdm;

  /// Vandermonde inverse V^{-1}: nodes → modes [N × N]
  RealMatrix m_vdmInv;

  /// V[i,0] value for the constant mode (1.0 for Legendre/quads, ~0.354 for Jacobi/triangles)
  CFreal m_vdm0;

  /// mode degree for each mode k (p_k in exp(-ζ*p_k²))
  std::vector<CFuint> m_modeDegrees;

  /// saved original nodal states [N × nbrEqs]
  RealMatrix m_origStates;

  /// modal coefficients û = V^{-1} * u [N × nbrEqs]
  RealMatrix m_modalCoeffs;

  /// temporary filtered nodal states [N × nbrEqs]
  RealMatrix m_filteredStates;

  /// temporary state vector for entropy computation
  RealVector m_tempState;

  /// filter weights per mode [N]
  RealVector m_filterWeights;

  /// filtered flux point states for bisection [maxNbrFlxPnts × nbrEqs]
  std::vector< RealVector > m_filteredFlxPntStates;

}; // class EntropyFilteringVS

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_EntropyFilteringVS_hh
