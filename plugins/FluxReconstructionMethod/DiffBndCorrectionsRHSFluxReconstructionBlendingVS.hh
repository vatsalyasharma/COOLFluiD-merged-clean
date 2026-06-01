// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionBlendingVS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionBlendingVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Scales diffusive boundary corrections by (1-alpha) for order blending.
 * Without this, boundary cells get full P1 diffusive BC contribution while
 * the interior diffusion is scaled by (1-alpha), creating a mismatch.
 *
 * @author Fix G7 — diffusive BC alpha-scaling for order blending
 */
class DiffBndCorrectionsRHSFluxReconstructionBlendingVS :
  public DiffBndCorrectionsRHSFluxReconstruction {

public:

  /// Constructor
  explicit DiffBndCorrectionsRHSFluxReconstructionBlendingVS(const std::string& name);

  /// Destructor
  virtual ~DiffBndCorrectionsRHSFluxReconstructionBlendingVS() {}

  /// Returns the DataSocket's that this command needs as sinks
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected:

  /// add the residual updates to the RHS, scaled by (1-alpha)
  virtual void updateRHS();

protected:

  /// socket for blending coefficient alpha
  Framework::DataSocketSink< CFreal > socket_alpha;

}; // class DiffBndCorrectionsRHSFluxReconstructionBlendingVS

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionBlendingVS_hh
