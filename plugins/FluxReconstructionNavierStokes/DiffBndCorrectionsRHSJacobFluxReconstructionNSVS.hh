#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionNSVS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionNSVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS for implicit schemes for NS.
   *
   * Adds BR2 damping and gradient variable transformation (setGradientVars)
   * to computeInterfaceFlxCorrection(). The existing computeBndGradTerms() already has
   * setGradientVars for the Jacobian perturbation path; this fix brings the same
   * transformation to the base RHS flux computation.
   *
   * @author Ray Vandenhoeck
   * @author Alexander Papen
   *
   */
class DiffBndCorrectionsRHSJacobFluxReconstructionNSVS : public DiffBndCorrectionsRHSJacobFluxReconstruction {

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSJacobFluxReconstructionNSVS(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSJacobFluxReconstructionNSVS();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * unset up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void unsetup();

protected: // functions

  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(CFreal& waveSpeedUpd);

  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);

  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);

  /**
   * compute the term for the gradient computation for the cell
   */
  virtual void computeCellGradTerm(RealMatrix& gradTerm);

  /**
   * compute the terms for the gradient computation for a face
   */
  virtual void computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR);

  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();

  /// compute the interface flux with BR2 damping and gradient variable transformation
  virtual void computeInterfaceFlxCorrection();

  /// set bnd face data — adds inverse characteristic length computation for BR2
  virtual void setBndFaceData(CFuint faceID);

protected: // data

    // vector for temporary storing the states of flx pnts (LEFT=cell, RIGHT=ghost)
    std::vector< std::vector< RealVector* > > m_tempStates;

    // vector for temporary storing the states of sol pnts
    std::vector< RealVector* > m_tempStatesSol;

    /// BR2 damping coefficient (from FluxReconstructionSolverData)
    CFreal m_dampCoeff;

    /// inverse characteristic lengths at flux points (computed per face)
    std::vector< CFreal > m_faceInvCharLengths;

    /// face flux point cell mapped coordinates (per local face index)
    Common::SafePtr< std::vector< std::vector< RealVector > > > m_faceFlxPntCellMappedCoords;

    /// temporary flux vector for BR2VS (average-of-fluxes) formulation
    RealVector m_tempFlux;

    /// whether to use BR2VS formulation
    bool m_useBR2VS;

}; // end of class DiffBndCorrectionsRHSJacobFluxReconstructionNSVS

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionNSVS_hh
