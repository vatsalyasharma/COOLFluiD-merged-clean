// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ExponentialModalFilterVS_hh
#define COOLFluiD_FluxReconstructionMethod_ExponentialModalFilterVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "MathTools/RealMatrix.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Exponential modal filter for shock capturing in Flux Reconstruction.
 *
 * Multiplies modal coefficients by sigma(eta) = exp(-alpha * eta^s)
 * where eta = k/N is the normalized mode number and s is the filter order.
 * This damps high-frequency modes while preserving the cell average (k=0).
 *
 * The filter matrix is precomputed in setup() as:
 *   FilterMatrix = V * diag(sigma_0, ..., sigma_{N-1}) * V^{-1}
 * and applied cell-locally to solution states each iteration.
 *
 * Placed in the LimiterCom slot (postProcessSolutionImpl), runs after
 * each time step update.
 *
 * Reference: GFR (NASA Glenn) exponential filter approach.
 */
class ExponentialModalFilterVS : public FluxReconstructionSolverCom {
public:

  /// Constructor
  explicit ExponentialModalFilterVS(const std::string& name);

  /// Destructor
  virtual ~ExponentialModalFilterVS();

  /// Defines the Config Option's of this class
  static void defineConfigOptions(Config::OptionList& options);

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Configures the command
  virtual void configure(Config::ConfigArgs& args);

  /// Execute the filter on all cells
  void execute();

protected: // data

  /// filter strength coefficient (alpha in exp(-alpha * eta^s))
  CFreal m_filterStrength;

  /// filter order (s in exp(-alpha * eta^s)), controls sharpness
  CFuint m_filterOrder;

  /// show rate for printing filter statistics
  CFuint m_showRate;

  /// whether the filter is active (DynamicOption for runtime control)
  bool m_active;

  /// precomputed filter matrix: V * diag(sigma) * V^{-1}
  RealMatrix m_filterMatrix;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// number of equations
  CFuint m_nbrEqs;

  /// number of solution points per cell
  CFuint m_nbrSolPnts;

  /// number of cells filtered this iteration
  CFuint m_nbFiltered;

  /// total number filtered (MPI-reduced)
  CFuint m_totalNbFiltered;

  /// temporary storage for original states (nbrSolPnts x nbrEqs)
  RealMatrix m_origStates;

}; // class ExponentialModalFilterVS

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ExponentialModalFilterVS_hh
