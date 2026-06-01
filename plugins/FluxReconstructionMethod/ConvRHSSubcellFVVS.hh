// Copyright (C) 2026
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvRHSSubcellFVVS_hh
#define COOLFluiD_FluxReconstructionMethod_ConvRHSSubcellFVVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective RHS with subcell FV blending (Hennemann-style).
 *
 * Both volume and surface terms are blended:
 *   R = (1-alpha)*R_FR + alpha*R_FV
 * where R_FV includes interior subcell faces and boundary subcell
 * faces (using the inter-element Riemann flux).
 */
class ConvRHSSubcellFVVS : public ConvRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvRHSSubcellFVVS(const std::string& name);

  /// Destructor
  virtual ~ConvRHSSubcellFVVS() {}

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
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvRHSSubcellFVVS_hh
