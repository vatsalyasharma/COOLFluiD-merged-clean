// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionNSBlendingVS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionNSBlendingVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffRHSFluxReconstructionBlendingVS.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// NS-specific diffusive RHS with order blending.
/// Inherits alpha scaling from DiffRHSFluxReconstructionBlendingVS
/// and adds NS-specific BR2 damping, viscosity-based wave speed,
/// and composition setting (same overrides as DiffRHSFluxReconstructionNS).
///
/// @author Rayan Dhib
class DiffRHSFluxReconstructionNSBlendingVS : public DiffRHSFluxReconstructionBlendingVS {

public: // functions

  /// Constructor
  explicit DiffRHSFluxReconstructionNSBlendingVS(const std::string& name);

  /// Destructor
  virtual ~DiffRHSFluxReconstructionNSBlendingVS() {}

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: //functions

  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);

  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();

  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();

 protected: //data

  /// matrix to store the state terms needed for the gradients for left neighbor
  RealMatrix m_tempGradTermL;

  /// matrix to store the state terms needed for the gradients for right neighbor
  RealMatrix m_tempGradTermR;

  /// diffusive variable set (NS-specific downcast)
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_diffusiveVarSet;

  /// element states of the left neighbor in the correct format
  std::vector< RealVector* > m_tempStatesL;

  /// element states of the right neighbor in the correct format
  std::vector< RealVector* > m_tempStatesR;

  /// BR2 damping coefficient
  CFreal m_dampCoeff;

  /// temporary flux vector for BR2VS (average-of-fluxes) formulation
  RealVector m_tempFlux;

  /// whether to use BR2VS formulation
  bool m_useBR2VS;

}; // class DiffRHSFluxReconstructionNSBlendingVS

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionNSBlendingVS_hh
