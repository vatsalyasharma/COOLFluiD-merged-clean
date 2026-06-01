#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQConsToRhoivtTvInRhoivtTvVS_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQConsToRhoivtTvInRhoivtTvVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "MathTools/MatrixInverter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative
 * to RhoivtTv variables, computed by building the forward Jacobian
 * dCons/dRhoivtTv and numerically inverting it.
 *
 * The forward matrix (dCons/dRhoivtTv) reuses the exact logic from
 * Euler2DNEQRhoivtTvToConsInRhoivtTv. The inverse gives the needed
 * dRhoivtTv/dCons transformation matrix for Newton's method.
 *
 * @author Vatsalya Sharma
 */
class Euler2DNEQConsToRhoivtTvInRhoivtTvVS : public Framework::VarSetMatrixTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;

  /**
   * Default constructor without arguments
   */
  Euler2DNEQConsToRhoivtTvInRhoivtTvVS(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQConsToRhoivtTvInRhoivtTvVS();

 /**
   * Set the transformation matrix from a given state
   */
  void setMatrix(const RealVector& state);

private:

  /**
   * Set the flag telling if the transformation is an identity one
   * @pre this method must be called during set up
   */
  bool getIsIdentityTransformation() const
  {
    return false;
  }

private: //data

  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;

  /// Vector storing the elemental composition
  RealVector _ys;

  /// array with all different vibrational dimensional temperatures
  RealVector _tvDim;

  /// array with all different vibrational dimensional energies
  RealVector _evDim;

  /// array to store density, enthalpy and energy
  RealVector _dhe;

  /// temporary matrix for the forward transformation (dCons/dRhoivtTv)
  RealMatrix _fwdMatrix;

  /// matrix inverter
  MathTools::MatrixInverter* _inverter;

}; // end of class Euler2DNEQConsToRhoivtTvInRhoivtTvVS

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQConsToRhoivtTvInRhoivtTvVS_hh
