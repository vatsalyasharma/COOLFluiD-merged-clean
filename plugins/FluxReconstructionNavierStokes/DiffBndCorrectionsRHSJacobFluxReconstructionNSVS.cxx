#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionNavierStokes/DiffBndCorrectionsRHSJacobFluxReconstructionNSVS.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionNSVS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
DiffBndCorrectionsRHSJacobNSVSFluxReconstructionProvider("DiffBndCorrectionsRHSJacobNSVS");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::DiffBndCorrectionsRHSJacobFluxReconstructionNSVS(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstruction(name),
  m_tempStates(),
  m_tempStatesSol(),
  m_dampCoeff(),
  m_faceInvCharLengths(),
  m_faceFlxPntCellMappedCoords(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::~DiffBndCorrectionsRHSJacobFluxReconstructionNSVS()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::setup()
{
  DiffBndCorrectionsRHSJacobFluxReconstruction::setup();

  m_tempStates.resize(2);
  m_tempStates[LEFT].resize(m_nbrFaceFlxPnts);
  m_tempStates[RIGHT].resize(m_nbrFaceFlxPnts);
  m_tempStatesSol.resize(m_nbrSolPnts);

  // BR2 damping coefficient
  m_dampCoeff = getMethodData().getDiffDampCoefficient();

  // get the local FR element data for cell mapped coordinates
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  m_faceFlxPntCellMappedCoords = frLocalData[0]->getFaceFlxPntCellMappedCoords();

  // resize inverse characteristic lengths
  m_faceInvCharLengths.resize(m_nbrFaceFlxPnts);

  // BR2VS config
  m_useBR2VS = getMethodData().getUseBR2VS();
  m_tempFlux.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::unsetup()
{
  DiffBndCorrectionsRHSJacobFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  CFreal visc = 1.0;
  /// @todo needs to be changed for non-NS
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  const CFreal dynVisc = navierStokesVarSet->getCurrDynViscosity();

  const CFreal factorPr = min(navierStokesVarSet->getModel().getPrandtl(),1.0);
  cf_assert(factorPr>0.0);

  waveSpeedUpd = 0.0;
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
    const CFreal rho = navierStokesVarSet->getDensity(*m_cellStatesFlxPnt[iFlx]);
    visc = dynVisc/rho/factorPr;

    // transform update states to physical data to calculate eigenvalues
    waveSpeedUpd += visc*jacobXJacobXIntCoef/m_cellVolume;
  }

}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::setBndFaceData(CFuint faceID)
{
  // call parent to set up all base geometric data
  DiffBndCorrectionsRHSJacobFluxReconstruction::setBndFaceData(faceID);

  // compute Jacobian determinants of the cell at flux point mapped coordinates
  std::valarray<CFreal> jacobDets =
    m_intCell->computeGeometricShapeFunctionJacobianDeterminant(
      (*m_faceFlxPntCellMappedCoords)[m_orient]);

  // compute inverse characteristic lengths for BR2 damping
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_faceInvCharLengths[iFlx] = m_faceJacobVecAbsSizeFlxPnts[iFlx] / jacobDets[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::computeInterfaceFlxCorrection()
{
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  // Collect state data pointers for gradient variable transformation
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx]  = m_cellStatesFlxPnt[iFlx]->getData();   // cell side
    m_tempStates[RIGHT][iFlx] = m_flxPntGhostSol[iFlx]->getData();     // ghost side
  }

  // Transform state variables to gradient variables (e.g., rho_i -> y_i for NEQ)
  // Reuses m_gradTermFace and m_ghostGradTerm from the Jacobian base class
  navierStokesVarSet->setGradientVars(m_tempStates[LEFT],  m_gradTermFace, m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(m_tempStates[RIGHT], m_ghostGradTerm, m_nbrFaceFlxPnts);

  // Loop over flux points to compute BR2 diffusive interface flux
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // compute the average solution and gradient (boundary: cell + ghost)
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[iFlxPnt][iVar])
                           + *(m_flxPntGhostGrads[iFlxPnt][iVar])) / 2.0;

      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[iFlxPnt]))[iVar]
                       + (*(m_flxPntGhostSol[iFlxPnt]))[iVar]) / 2.0;
    }

    // BR2 damping: avgGrad -= dampFactor * (y_cell - y_ghost) * n
    const CFreal dampFactor = m_dampCoeff * m_faceInvCharLengths[iFlxPnt];

    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      const RealVector dGradVarXNormal = (m_gradTermFace(iGrad, iFlxPnt)
                                         - m_ghostGradTerm(iGrad, iFlxPnt))
                                        * m_unitNormalFlxPnts[iFlxPnt];
      *m_avgGrad[iGrad] -= dampFactor * dGradVarXNormal;
    }

    if (m_useBR2VS)
    {
      // BR2VS: average-of-fluxes formulation
      // Cell-side flux: F_v(u_cell, avgGrad_BR2)
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_avgSol[iVar] = (*(m_cellStatesFlxPnt[iFlxPnt]))[iVar];
      }
      prepareFluxComputation();
      computeFlux(m_avgSol, m_avgGrad, m_unitNormalFlxPnts[iFlxPnt], 0, m_tempFlux);

      // Ghost-side flux: F_v(u_ghost, avgGrad_BR2)
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_avgSol[iVar] = (*(m_flxPntGhostSol[iFlxPnt]))[iVar];
      }
      prepareFluxComputation();
      computeFlux(m_avgSol, m_avgGrad, m_unitNormalFlxPnts[iFlxPnt], 0, m_flxPntRiemannFlux[iFlxPnt]);

      // Average of the two fluxes
      m_flxPntRiemannFlux[iFlxPnt] = 0.5 * (m_tempFlux + m_flxPntRiemannFlux[iFlxPnt]);
    }
    else
    {
      // Standard BR2: F_v(avg_u, avgGrad_BR2)
      prepareFluxComputation();
      computeFlux(m_avgSol, m_avgGrad, m_unitNormalFlxPnts[iFlxPnt], 0, m_flxPntRiemannFlux[iFlxPnt]);
    }

    // scale by face Jacobian (boundary: single cell, scalar)
    m_cellFlx[iFlxPnt] = m_flxPntRiemannFlux[iFlxPnt] * m_faceJacobVecSizeFlxPnts[iFlxPnt];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_cellStatesFlxPnt[iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_flxPntGhostSol[iFlx]->getData());
  }

  navierStokesVarSet->setGradientVars(m_tempStates[LEFT],gradTerm,m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(m_tempStates[RIGHT],ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_cellStatesFlxPnt2[iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_flxPntGhostSol[iFlx]->getData());
  }

  navierStokesVarSet->setGradientVars(m_tempStates[LEFT],gradTerm,m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(m_tempStates[RIGHT],ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::computeCellGradTerm(RealMatrix& gradTerm)
{
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStatesSol[iSol] = ((*m_cellStates)[iSol]->getData());
  }

  navierStokesVarSet->setGradientVars(m_tempStatesSol,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_pertCellStatesFlxPnt[LEFT][iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_pertCellStatesFlxPnt[RIGHT][iFlx]->getData());
  }

  navierStokesVarSet->setGradientVars(m_tempStates[LEFT],gradTermL,m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(m_tempStates[RIGHT],gradTermR,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSVS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
  const bool effectivePerturb = this->getMethodData().freezeDiffCoeff() ? false : isPerturb;
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  navierStokesVarSet->setComposition(m_avgSol, effectivePerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
