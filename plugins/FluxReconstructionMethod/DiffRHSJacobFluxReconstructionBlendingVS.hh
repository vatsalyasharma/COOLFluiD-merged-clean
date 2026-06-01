// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionBlendingVS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionBlendingVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffRHSJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// Diffusive RHS + Jacobian with order blending: scales both the diffusive
/// residual and its Jacobian by (1-alpha) per cell, consistent with
/// convective order blending (ConvRHSJacobFluxReconstructionBlending).
///
/// Alpha is frozen during FD perturbation, so d[(1-a)*R]/du = (1-a)*dR/du.
///
/// @author Vatsalya Sharma
class DiffRHSJacobFluxReconstructionBlendingVS : public DiffRHSJacobFluxReconstruction {

public:

  explicit DiffRHSJacobFluxReconstructionBlendingVS(const std::string& name);

  virtual ~DiffRHSJacobFluxReconstructionBlendingVS() {}

  virtual void execute();

  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected:

  virtual void computeBothJacobsDiffFaceTerm();

  virtual void computeOneJacobDiffFaceTerm(const CFuint side);

protected:

  /// socket for blending coefficient alpha per cell
  Framework::DataSocketSink< CFreal > socket_alpha;

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionBlendingVS_hh
