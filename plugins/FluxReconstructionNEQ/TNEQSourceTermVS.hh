#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_TNEQSourceTermVS_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_TNEQSourceTermVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNEQ/CNEQSourceTermVS.hh"

#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * A command for adding the source term for TNEQ
 *
 * @author Ray Vandenhoeck
 * @author Firas Ben Ameur
 *
 */
class TNEQSourceTermVS : public CNEQSourceTermVS {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit TNEQSourceTermVS(const std::string& name);

  /**
   * Destructor.
   */
  ~TNEQSourceTermVS();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected:

  /**
   * get data required for source term computation
   */
  void getSourceTermData();

  /**
   * add the source term
   */
  void addSourceTerm(RealVector& resUpdates);

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set the adimensional vibrational temperatures
   */
  virtual void setVibTemperature(const RealVector& pdata,
				 const Framework::State& state,
				 RealVector& tvib);

  /// Compute non conservative term \f$ p_e \nabla \cdot v \f$.
  void computePeDivV(Framework::GeometricEntity *const element,
		     RealVector& source, RealMatrix& jacobian);

  /// Compute the energy transfer terms
  void computeSourceVT(RealVector& omegaTv, CFreal& omegaRad);

  /// Override: analytical source Jacobian using PLATO get_source_jac
  void getSToStateJacobian(const CFuint iState); //Vatsalya: same as 1T solution but for TCNEQ; numerical FD Jacobian could not capture exponential T-sensitivity of V-T relaxation (tau_VT~1e-8s), limiting CFL to 0.01

protected: // data

  /// PLATO Jacobian matrix (nEqs x nEqs in PLATO variable space)
  RealMatrix m_platoJacob;

  /// radiation source term
  CFdouble m_omegaRad;

  /// velocity divergence
  CFreal m_divV;

  /// electron pressure
  CFreal m_pe;

  /// array to store the energy exchange terms
  RealVector m_omegaTv;

  /// reference data
  Common::SafePtr<RealVector> m_refData;

}; // class TNEQSourceTermVS

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_TNEQSourceTermVS_hh
