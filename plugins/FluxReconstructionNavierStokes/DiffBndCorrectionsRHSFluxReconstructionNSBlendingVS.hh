#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNSBlendingVS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNSBlendingVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/DiffBndCorrectionsRHSFluxReconstructionNS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * NS boundary diffusive corrections with BR2 damping + order blending (1-alpha) scaling.
 * Inherits all BR2/setGradientVars from the NS parent; adds alpha-scaling in updateRHS().
 *
 * Provider: "DiffBndCorrectionsRHSNSBlendingVS"
 */
class DiffBndCorrectionsRHSFluxReconstructionNSBlendingVS :
  public DiffBndCorrectionsRHSFluxReconstructionNS {

public:

  /// Constructor
  explicit DiffBndCorrectionsRHSFluxReconstructionNSBlendingVS(const std::string& name);

  /// Destructor
  virtual ~DiffBndCorrectionsRHSFluxReconstructionNSBlendingVS() {}

  /// Returns the DataSocket's that this command needs as sinks
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected:

  /// add the residual updates to the RHS, scaled by (1-alpha) for order blending consistency
  virtual void updateRHS();

protected:

  /// socket for blending coefficient alpha
  Framework::DataSocketSink< CFreal > socket_alpha;

}; // class DiffBndCorrectionsRHSFluxReconstructionNSBlendingVS

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNSBlendingVS_hh
