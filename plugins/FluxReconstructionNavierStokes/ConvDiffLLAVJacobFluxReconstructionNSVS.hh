// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstructionNSVS_hh
#define COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstructionNSVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/ConvDiffLLAVJacobFluxReconstructionNSBaseVS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * LLAV-only Jacobian preconditioner for implicit FR.
 *
 * During frozen-Jacobian iterations (KeepLLAVInJacob=true), this class
 * replaces the stale full Jacobian with a freshly computed matrix containing
 * ONLY the LLAV diffusion Jacobian + M/dt diagonal.  This avoids the crash
 * caused by FD perturbation of convective fluxes at shocks, while keeping
 * M/dt current with the CFL ramp.
 *
 * When KeepLLAVInJacob=false (default) or the Jacobian is not frozen, the
 * class follows the FR-derived ConvDiffLLAVJacobFluxReconstructionNSBaseVS
 * path.
 *
 * @author Vatsalya Sharma
 */
class ConvDiffLLAVJacobFluxReconstructionNSVS : public ConvDiffLLAVJacobFluxReconstructionNSBaseVS {

public: // functions

  /// Constructor
  explicit ConvDiffLLAVJacobFluxReconstructionNSVS(const std::string& name);

  /// Destructor
  virtual ~ConvDiffLLAVJacobFluxReconstructionNSVS() {}

  /// Define config options
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure
  virtual void configure(Config::ConfigArgs& args);

  /// Set up private data
  virtual void setup();

  /// Execute: base RHS+Jacob, then optional LLAV-only second pass
  virtual void execute();

protected: // functions

  /// Second face loop: compute LLAV Jacobian + M/dt for both sides
  void computeLLAVOnlyBothJacobs();

  /// Second face loop: compute LLAV Jacobian + M/dt for one side (MPI boundary)
  void computeLLAVOnlyOneJacob(const CFuint side);

  /// Add M/dt diagonal block for one cell side
  void addTimeDiagonal(const CFuint side);

protected: // data

  /// DynamicOption: when true AND Jacobian is frozen, recompute LLAV+M/dt matrix
  bool m_keepLLAVInJacob;

  /// Diagonal values for M/dt (per solution point)
  std::vector< CFreal > m_diagValuesVS;

  /// Solution-point local coordinates (for Jacobian determinant computation)
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoordsVS;

}; // class

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstructionNSVS_hh
