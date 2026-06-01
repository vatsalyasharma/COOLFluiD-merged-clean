// Copyright (C) 2026
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvRHSJacobSubcellFVVS_hh
#define COOLFluiD_FluxReconstructionMethod_ConvRHSJacobSubcellFVVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective RHS+Jacobian with subcell FV blending (Hennemann-style).
 *
 * Both volume and surface terms are blended:
 *   R = (1-alpha)*R_FR + alpha*R_FV
 * where R_FV includes interior subcell faces and boundary subcell
 * faces (using the inter-element Riemann flux).
 */
class ConvRHSJacobSubcellFVVS : public ConvRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvRHSJacobSubcellFVVS(const std::string& name);

  /// Destructor
  virtual ~ConvRHSJacobSubcellFVVS() {}

  /// Execute processing actions
  void execute();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected: // functions

  /// add the updates to the wave speed (blended with alpha)
  void updateWaveSpeed();

  /// compute the divergence of the discontinuous flx (-divFD+divhFD) with subcell FV blending
  void computeDivDiscontFlx(std::vector< RealVector >& residuals);

  /// compute subcell FV residual for current cell (volume term only)
  void computeSubcellFVResidual();

  /// compute the Jacobian contribution for the volume term with subcell FV blending
  void computeJacobConvCorrection();

  /// compute analytical subcell FV Jacobian and add to accumulator
  void computeAnalyticalSubcellFVJacobianVS(Framework::BlockAccumulator& acc,
                                           CFreal resFactor, CFreal alpha);

  /// compute physical data Jacobian ∂pdata/∂U by FD of computePhysicalData
  void computePhysicalDataJacobianVS(Framework::State& state,
                                    const RealVector& pdataBase,
                                    RealMatrix& dPdU);

  /// compute analytical AUSM+ flux Jacobian ∂F/∂U for a subcell interface.
  /// Fills dFluxdUL and dFluxdUR (nEqs × nEqs each).
  void computeAUSMPlusFluxJacobianVS(
    const RealVector& pdataL, const RealVector& pdataR,
    const RealVector& unitNormal,
    const RealMatrix& dPdUL, const RealMatrix& dPdUR,
    RealMatrix& dFluxdUL, RealMatrix& dFluxdUR);

protected: // data

  /// socket for output of the filtering (alpha)
  Framework::DataSocketSink< CFreal > socket_alpha;

  /// 1D solution point coordinates (reference)
  std::vector< CFreal > m_solPnts1D;

  /// 1D subcell interfaces (reference)
  std::vector< CFreal > m_subcellInterfaces1D;

  /// inverse 1D subcell widths
  std::vector< CFreal > m_invSubcellWidth1D;

  /// number of solution points per direction
  CFuint m_nbrSolPnts1D;

  /// whether element is tensor-product (quad/hex)
  bool m_tensorBased;

  /// subcell interface coordinates per direction (reference)
  std::vector< std::vector< RealVector > > m_subcellIfaceCoords;

  /// dimension list per direction (same size as coords)
  std::vector< std::vector< CFuint > > m_subcellIfaceDimList;

  /// left/right solution indices per interface
  std::vector< std::vector< CFuint > > m_subcellIfaceLeftIdx;
  std::vector< std::vector< CFuint > > m_subcellIfaceRightIdx;

  /// inverse widths for left/right subcells per interface
  std::vector< std::vector< CFreal > > m_subcellIfaceInvWidthL;
  std::vector< std::vector< CFreal > > m_subcellIfaceInvWidthR;

  /// subcell FV residual (per solution point)
  std::vector< RealVector > m_subcellFVRes;

  /// unit normal for subcell FV Riemann flux (reused across interfaces)
  RealVector m_subcellUnitNormal;

  /// threshold to skip FV computation when alpha is tiny
  CFreal m_alphaEps;

  /// boundary subcell solution point index for each flux point (tensor-product only)
  std::vector< CFuint > m_flxPntBndSolIdx;

  /// inverse subcell width for boundary subcell for each flux point
  std::vector< CFreal > m_flxPntBndInvWidth;

  // --- Analytical subcell FV Jacobian storage (AJ-1) ---

  /// physical data per solution point (cached for Jacobian)
  std::vector< RealVector > m_pdataPerSol;

  /// physical data Jacobian per solution point: ∂pdata/∂U [pdataSize × nEqs]
  std::vector< RealMatrix > m_dPdUPerSol;

  /// temporary perturbed physical data
  RealVector m_pdataPerturbed;

  /// physical data size
  CFuint m_pdataSize;

  /// flux Jacobian w.r.t. left/right conservative state [nEqs × nEqs]
  RealMatrix m_dFluxdUL;
  RealMatrix m_dFluxdUR;

  /// number of species equations (nEqs - dim - 1)
  CFuint m_nbSpecies;

  /// index where Euler equations start in the equation ordering
  CFuint m_eulerStartIdx;

  /// index of first scalar variable in physical data
  CFuint m_firstScalarPdata;

  /// FD epsilon for physical data Jacobian
  static constexpr CFreal m_pdataFDEps = 1.0e-7;
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvRHSJacobSubcellFVVS_hh
