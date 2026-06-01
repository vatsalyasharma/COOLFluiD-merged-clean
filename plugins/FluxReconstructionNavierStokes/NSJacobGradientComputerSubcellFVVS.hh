// Copyright (C) 2026
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_NSJacobGradientComputerSubcellFVVS_hh
#define COOLFluiD_FluxReconstructionMethod_NSJacobGradientComputerSubcellFVVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSJacobSubcellFVVS.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Daughterclass of ConvRHSJacobSubcellFVVS, needed to calculate the gradients
 * for implicit NS schemes with subcell FV blending (adds setGradientVars transformation).
 */
class NSJacobGradientComputerSubcellFVVS : public ConvRHSJacobSubcellFVVS {

public: // functions

  /// Constructor
  explicit NSJacobGradientComputerSubcellFVVS(const std::string& name);

  /// Destructor
  virtual ~NSJacobGradientComputerSubcellFVVS() {}

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: //functions

  /**
   * Compute the discontinuous contribution to the corrected gradients
   */
  virtual void computeGradients();

  /**
   * Compute the correction part of the corrected gradient
   */
  virtual void computeGradientFaceCorrections();

protected: //data

  /// diffusive variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_diffusiveVarSet;

  /// matrix to store the state terms needed for the gradients (y_i, u, v, T) inside element
  RealMatrix m_tempGradTerm;

  /// matrix to store the state terms needed for the gradients (y_i, u, v, T) for left neighbor
  RealMatrix m_tempGradTermL;

  /// matrix to store the state terms needed for the gradients (y_i, u, v, T) for right neighbor
  RealMatrix m_tempGradTermR;

  /// element states within an element in the correct format
  std::vector< RealVector* > m_tempStates;

  /// element states of the left neighbor in the correct format
  std::vector< RealVector* > m_tempStatesL;

  /// element states of the right neighbor in the correct format
  std::vector< RealVector* > m_tempStatesR;

}; // class NSJacobGradientComputerSubcellFVVS

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NSJacobGradientComputerSubcellFVVS_hh
