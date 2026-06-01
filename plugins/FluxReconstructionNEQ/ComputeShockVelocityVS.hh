#ifndef COOLFluiD_FluxReconstructionNEQ_ComputeShockVelocityVS_hh
#define COOLFluiD_FluxReconstructionNEQ_ComputeShockVelocityVS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Computes the Rankine-Hugoniot shock velocity diagnostic at boundary faces.
 *
 * For each face on the configured TRS (ShockDown), extrapolates the interior
 * solution to the face flux points (using the FR polynomial), then computes:
 *   Vs = (rho_fs * u_n_fs - rho_int * u_n_int) / (rho_fs - rho_int)
 * at each flux point, and averages over the face.
 *
 * Vs = 0 means the shock is at steady state at this position.
 * Vs > 0 means the shock wants to move toward the freestream.
 * Vs < 0 means the shock wants to move toward the wall.
 *
 * Stores Vs in a per-state socket for VTU output via DataHandleOutput.
 *
 * @author COOLFluiD Shock Fitting Extension
 */
class ComputeShockVelocityVS : public Framework::DataProcessingCom {

public:

  static void defineConfigOptions(Config::OptionList& options);

  ComputeShockVelocityVS(const std::string& name);
  ~ComputeShockVelocityVS();

  static std::string getClassName() { return "ComputeShockVelocityVS"; }

  void setup();
  void unsetup();
  void configure(Config::ConfigArgs& args);

  std::vector<Common::SafePtr<Framework::BaseDataSocketSource>> providesSockets();
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink>> needsSockets();

protected:

  void executeOnTrs();

protected:

  /// Socket: shock velocity per state (0 for non-boundary cells)
  Framework::DataSocketSource<CFreal> socket_shockVelocity;

  /// Socket: mass flux imbalance per state (0 for non-boundary cells)
  Framework::DataSocketSource<CFreal> socket_massFluxImbalance;

  /// Sink: states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// Sink: face Jacobian vector sizes at face flux points
  Framework::DataSocketSink<std::vector<CFreal>> socket_faceJacobVecSizeFaceFlxPnts;

  /// FR solver data
  Common::SafePtr<FluxReconstructionSolverData> m_frData;

  /// Face builder
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder>> m_faceBuilder;

  /// Cell builder
  Common::SafePtr<Framework::GeometricEntityPool<FluxReconstructionMethod::CellToFaceGEBuilder>> m_cellBuilder;

  /// Freestream state values
  std::vector<CFreal> m_fsDef;

  /// Freestream total density, velocity
  CFreal m_rhoFs, m_uFs, m_vFs;

  /// Number of solution points per cell
  CFuint m_nbrSolPnts;

  /// Number of face flux points
  CFuint m_nbrFaceFlxPnts;

  /// Number of species
  CFuint m_nbSpecies;

  /// Polynomial values at flux points (for extrapolation from sol points)
  Common::SafePtr<std::vector<std::vector<CFreal>>> m_solPolyValsAtFlxPnts;

  /// Face flux point connectivity (orient → flux point → element flux point index)
  Common::SafePtr<std::vector<std::vector<CFuint>>> m_faceFlxPntConn;

  /// Face mapped coordinate direction
  Common::SafePtr<std::vector<CFint>> m_faceMappedCoordDir;

  /// Flux point local coordinates on face
  Common::SafePtr<std::vector<RealVector>> m_flxLocalCoords;

  /// Extrapolated state at each face flux point
  std::vector<RealVector> m_cellStatesFlxPnt;

  /// Unit normal at each face flux point
  std::vector<RealVector> m_unitNormalFlxPnts;

  /// Current face orientation
  CFuint m_orient;

}; // class ComputeShockVelocityVS

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
