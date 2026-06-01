#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNSVS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNSVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSFluxReconstruction.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS for NS.
   *
   * Adds BR2 damping and gradient variable transformation (setGradientVars)
   * that the base class lacks. Interior faces have these in DiffRHSFluxReconstructionNSVS;
   * this class provides the boundary-face equivalent.
   *
   * @author Ray Vandenhoeck
   * @author Alexander Papen
   *
   */
class DiffBndCorrectionsRHSFluxReconstructionNSVS : public DiffBndCorrectionsRHSFluxReconstruction {

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSFluxReconstructionNSVS(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSFluxReconstructionNSVS();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: // functions

  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(CFreal& waveSpeedUpd);

  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();

  /// compute the interface flux with BR2 damping and gradient variable transformation
  virtual void computeInterfaceFlxCorrection();

  /// set bnd face data — adds inverse characteristic length computation for BR2
  virtual void setBndFaceData(CFuint faceID);

protected: //data

  /// matrix to store the gradient variables for cell-side flux point states
  RealMatrix m_tempGradTerm;

  /// matrix to store the gradient variables for ghost-side flux point states
  RealMatrix m_tempGradTermGhost;

  /// diffusive variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_diffusiveVarSet;

  /// cell-side states at flux points in the correct format for setGradientVars
  std::vector< RealVector* > m_tempStates;

  /// ghost-side states at flux points in the correct format for setGradientVars
  std::vector< RealVector* > m_tempStatesGhost;

  /// BR2 damping coefficient (from FluxReconstructionSolverData)
  CFreal m_dampCoeff;

  /// inverse characteristic lengths at flux points (computed per face in setBndFaceData)
  std::vector< CFreal > m_faceInvCharLengths;

  /// face flux point cell mapped coordinates (per local face index, for Jacobian determinant)
  Common::SafePtr< std::vector< std::vector< RealVector > > > m_faceFlxPntCellMappedCoords;

  /// temporary flux vector for BR2VS (average-of-fluxes) formulation
  RealVector m_tempFlux;

  /// whether to use BR2VS formulation
  bool m_useBR2VS;

}; // end of class DiffBndCorrectionsRHSFluxReconstructionNSVS

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNSVS_hh
