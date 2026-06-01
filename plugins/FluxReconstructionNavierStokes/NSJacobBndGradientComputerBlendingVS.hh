// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_NSJacobBndGradientComputerBlendingVS_hh
#define COOLFluiD_FluxReconstructionMethod_NSJacobBndGradientComputerBlendingVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSJacobFluxReconstructionBlending.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Daughterclass of ConvBndCorrectionsRHSJacobFluxReconstructionBlending,
 * needed to calculate the bnd gradients for implicit schemes with order
 * blending for NS.
 *
 * Combines the NS diffusive varset gradient logic from NSJacobBndGradientComputer
 * with the blending+Jacobian machinery from ConvBndCorrectionsRHSJacobFluxReconstructionBlending.
 *
 * @author Vatsalya Sharma
 */
class NSJacobBndGradientComputerBlendingVS : public ConvBndCorrectionsRHSJacobFluxReconstructionBlending {

public: // functions

  /// Constructor
  explicit NSJacobBndGradientComputerBlendingVS(const std::string& name);

  /// Destructor
  virtual ~NSJacobBndGradientComputerBlendingVS() {}

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: //functions

  /**
   * Compute the correction part of the corrected gradient for the bnd face
   */
  virtual void computeGradientBndFaceCorrections();

protected: //data

  /// diffusive variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_diffusiveVarSet;

  /// Vector transformer from update to solution variables
  Common::SafePtr<Framework::VarSetTransformer> m_updateToSolutionVecTrans;

  /// matrix to store the state terms needed for the gradients (p, u, v, T)
  RealMatrix m_tempGradTerm;

  /// matrix to store the ghost state terms needed for the gradients (p, u, v, T)
  RealMatrix m_tempGradTermGhost;

  /// element states within an element in the correct format
  std::vector< RealVector* > m_tempStates;

  /// element ghost states in the correct format
  std::vector< RealVector* > m_tempStatesGhost;

}; // class NSJacobBndGradientComputerBlendingVS

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NSJacobBndGradientComputerBlendingVS_hh
