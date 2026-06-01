// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionBlendingVS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionBlendingVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// Diffusive RHS with order blending: scales diffusive contributions
/// by (1-alpha) per cell, consistent with convective order blending.
///
/// When alpha=0 (smooth region), full P1 diffusion is applied.
/// When alpha=1 (shock region), diffusion is suppressed (P0 convection
/// provides sufficient numerical dissipation via the Riemann solver).
///
/// Mathematically: the P0 volume divergence of diffusive flux is zero
/// (constant polynomial => zero derivative), so the blended diffusive
/// residual is (1-alpha)*R_diff_P1 + alpha*0 = (1-alpha)*R_diff_P1.
///
/// @author Rayan Dhib
class DiffRHSFluxReconstructionBlendingVS : public DiffRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit DiffRHSFluxReconstructionBlendingVS(const std::string& name);

  /// Destructor
  virtual ~DiffRHSFluxReconstructionBlendingVS() {}

  /// Execute processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected: //data

  /// socket for blending coefficient alpha per cell
  Framework::DataSocketSink< CFreal > socket_alpha;

}; // class DiffRHSFluxReconstructionBlendingVS

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionBlendingVS_hh
