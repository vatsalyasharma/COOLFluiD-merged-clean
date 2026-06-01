#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionPrim_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHD3DProjectionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 3D for projection scheme
 * for primitive variables
 *
 * @author Andrea Lani
 * @author Radka Keslerova
 */

class MHD3DProjectionPrim : public MHD3DProjectionVarSet {
public: //function

  /**
   * Constructor
   * @see MHD3DProjection
   */
  MHD3DProjectionPrim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~MHD3DProjectionPrim();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();

  /**
   * Get extra variable names
   */
  virtual std::vector<std::string> getExtraVarNames() const;

  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const;

  /**
   * Set the jacobians
   */
  virtual void computeJacobians();

  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
			     RealMatrix& jacobMin,
			     RealVector& eValues,
			     const RealVector& normal);
  
  /**
   * Set the matrix of the right and left eigenvectors
   * and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
					 RealMatrix& leftEv,
					 RealVector& eValues,
					 const RealVector& normal);
  
  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  virtual void computePhysicalData(const Framework::State& state,
				   RealVector& data);
  
  /**
   * Set the total magnetic field and energy values
   */
  virtual void setDimensionalValuesPlusExtraValues(const Framework::State& state,
						   RealVector& result,
						   RealVector& extra);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  virtual void computeStateFromPhysicalData(const RealVector& data,
					    Framework::State& state);
  
  
private: // helper function

  /**
   * Set the rotation matrices
   */
  void setRotationMatrices(const RealVector& normals);

private: // data
  // Rotation matrix
  RealMatrix _rm;

  // Inverse rotation matrix
  RealMatrix _rmInv;

  // Temporary matrix of right eigenvectors
  RealMatrix _rightEv;

  // Temporary matrix of left eigenvectors
  RealMatrix _leftEv;

}; // end of class MHD3DProjectionPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif //COOLFluiD_Physics_MHD_MHD3DProjectionPrim_hh
