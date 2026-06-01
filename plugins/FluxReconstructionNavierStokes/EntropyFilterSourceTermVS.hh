#ifndef COOLFluiD_FluxReconstructionMethod_EntropyFilterSourceTermVS_hh
#define COOLFluiD_FluxReconstructionMethod_EntropyFilterSourceTermVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/StdSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Implicit modal dissipation source term for FR.
 *
 * Adds S_filter(u) = -kappa * (I - F) * u to the residual, where
 * F = V * diag(sigma_k) * V^{-1} is the modal filter matrix with
 * sigma_k = exp(-zeta * p_k^2).
 *
 * The filter strength zeta is computed self-consistently from the
 * modal energy ratio of the current cell states:
 *   eta = E_high / (E_high + E_mean)
 *   zeta = zetaMax * clamp((eta - threshold) / (saturation - threshold))
 *
 * This makes modal dissipation visible to Newton's FD Jacobian.
 * At steady state: R_PDE - kappa*(I-F)*u = 0. Smooth cells (zeta=0)
 * get the exact PDE solution. Shock cells get damped high-frequency content.
 *
 * @author Vatsalya Sharma
 */
class EntropyFilterSourceTermVS : public StdSourceTerm {
public:

  static void defineConfigOptions(Config::OptionList& options);

  explicit EntropyFilterSourceTermVS(const std::string& name);

  virtual ~EntropyFilterSourceTermVS();

  virtual void setup();

  virtual void unsetup();

  virtual void addSourceTerm(RealVector& resUpdates);

  virtual void getSourceTermData();

protected: // data

  /// Vandermonde matrix [N x N]
  RealMatrix m_vdm;

  /// Vandermonde inverse [N x N]
  RealMatrix m_vdmInv;

  /// mode polynomial degree for each mode k
  std::vector<CFuint> m_modeDegrees;

  /// work array: modal coefficients [N x nbrEqs]
  RealMatrix m_modalCoeffs;

  /// work array: dissipation weights d_k = 1 - exp(-zeta * p_k^2) [N]
  RealVector m_dissipWeights;

  /// filter relaxation rate (DynamicOption)
  CFreal m_kappa;

  /// maximum filter strength (exp(-zetaMax) ~ 0 for P1 modes)
  CFreal m_zetaMax;

  /// modal energy threshold below which no filtering is applied
  CFreal m_threshold;

  /// modal energy saturation level (full filtering applied)
  CFreal m_saturation;

  /// enable/disable at runtime (DynamicOption)
  bool m_active;

  /// log every N iterations
  CFuint m_showRate;

  /// number of dimensions
  CFuint m_nbDims;

  /// iteration counter for active cell logging
  mutable CFuint m_nActiveCells;

}; // class EntropyFilterSourceTermVS

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_EntropyFilterSourceTermVS_hh
