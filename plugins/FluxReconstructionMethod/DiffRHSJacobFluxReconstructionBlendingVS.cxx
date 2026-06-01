// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Diffusive Jacobian class for order blending.
// SFV-6: Per Rueda-Ramirez et al. (2021, resistive MHD), diffusion is NOT blended:
//   "the visco-resistive terms are discretized only with the high-order DGSEM method"
// All (1-alpha) scaling removed. Full FR diffusion kept in all cells regardless of alpha.
// Alpha socket retained for provider routing compatibility.

#include "Common/PE.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/DiffRHSJacobFluxReconstructionBlendingVS.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffRHSJacobFluxReconstructionBlendingVS,
                       FluxReconstructionSolverData,
                       FluxReconstructionModule >
diffRHSJacobBlendingVSFluxReconstructionProvider("DiffRHSJacobBlendingVS");

//////////////////////////////////////////////////////////////////////////////

DiffRHSJacobFluxReconstructionBlendingVS::DiffRHSJacobFluxReconstructionBlendingVS(const std::string& name) :
  DiffRHSJacobFluxReconstruction(name),
  socket_alpha("alpha")
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
DiffRHSJacobFluxReconstructionBlendingVS::needsSockets()
{
  // Get parent sockets
  std::vector< Common::SafePtr< BaseDataSocketSink > > result =
    DiffRHSJacobFluxReconstruction::needsSockets();
  // Add alpha socket for blending coefficient
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionBlendingVS::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "DiffRHSJacobFluxReconstructionBlendingVS::execute()\n");

  // Get blending coefficient alpha
  DataHandle< CFreal > alphaData = socket_alpha.getDataHandle();

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;

  // reset cell flags
  for (CFuint iCell = 0; iCell < m_cellFlags.size(); ++iCell)
  {
    m_cellFlags[iCell] = false;
  }

  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;

  // get the geodata of the cell builders and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCBL = m_cellBuilders[LEFT]->getDataGE();
  geoDataCBL.trs = cells;
  CellToFaceGEBuilder::GeoData& geoDataCBR = m_cellBuilders[RIGHT]->getDataGE();
  geoDataCBR.trs = cells;

  //// Loop over faces to calculate fluxes and interface fluxes in the flux points

  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    CFLog(VERBOSE, "Orient = " << m_orient << "\n");
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
    m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();

      // compute volume
      m_cellVolume[LEFT] = m_cells[LEFT]->computeVolume();
      m_cellVolume[RIGHT] = m_cells[RIGHT]->computeVolume();

      cf_assert(m_cellVolume[LEFT] > 0.0);
      cf_assert(m_cellVolume[RIGHT] > 0.0);

      // if one of the neighbouring cells is parallel updatable, compute the correction flux
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
        // build the neighbouring cells
        const CFuint cellIDL = m_face->getNeighborGeo(LEFT)->getID();
        geoDataCBL.idx = cellIDL;
        m_cells[LEFT] = m_cellBuilders[LEFT ]->buildGE();
        const CFuint cellIDR = m_face->getNeighborGeo(RIGHT)->getID();
        geoDataCBR.idx = cellIDR;
        m_cells[RIGHT] = m_cellBuilders[RIGHT]->buildGE();

        // set the face data
        setFaceData(m_face->getID());

        // compute the left and right states and gradients in the flx pnts
        computeFlxPntStatesAndGrads();

        // compute FI
        computeInterfaceFlxCorrection();

        // compute the wave speed updates
        computeWaveSpeedUpdates(m_waveSpeedUpd);

        // update the wave speed
        updateWaveSpeed();

        // === LEFT face correction (SFV-6: no alpha scaling on diffusion) ===
        computeCorrection(LEFT, m_divContFlxL);
        m_divContFlx = m_divContFlxL;

        // update RHS
        updateRHS();

        // === RIGHT face correction (SFV-6: no alpha scaling on diffusion) ===
        computeCorrection(RIGHT, m_divContFlxR);
        m_divContFlx = m_divContFlxR;

        // update RHS
        updateRHS();

        // === Cell volume LEFT (SFV-6: no alpha scaling on diffusion) ===
        if (!m_cellFlags[cellIDL] && (*m_states[LEFT ])[0]->isParUpdatable())
        {
          computeUnpertCellDiffResiduals(LEFT);
          m_unpertAllCellDiffRes[cellIDL] = m_unpertCellDiffRes[LEFT];

          // update RHS
          updateRHSUnpertCell(LEFT);
        }

        // === Cell volume RIGHT (SFV-6: no alpha scaling on diffusion) ===
        if (!m_cellFlags[cellIDR] && (*m_states[RIGHT])[0]->isParUpdatable())
        {
          computeUnpertCellDiffResiduals(RIGHT);
          m_unpertAllCellDiffRes[cellIDR] = m_unpertCellDiffRes[RIGHT];

          // update RHS
          updateRHSUnpertCell(RIGHT);
        }

        // === Jacobian (gated by FreezeJacob) ===
        {
          const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
          const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
          const CFuint interval = iter - iterFreeze;
          if (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0)
          {
            // get all the faces neighbouring the cells
            m_faces[LEFT ] = m_cells[LEFT ]->getNeighborGeos();
            m_faces[RIGHT] = m_cells[RIGHT]->getNeighborGeos();

            // set the local indexes of the other faces than the current faces
            setOtherFacesLocalIdxs();

            // make a back up of the grads and put the perturbed and unperturbed corrections in the correct format
            for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
            {
              for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
              {
                m_cellGradsBackUp[LEFT][iState][iVar] = (*m_cellGrads[LEFT][iState])[iVar];
                m_cellGradsBackUp[RIGHT][iState][iVar] = (*m_cellGrads[RIGHT][iState])[iVar];

                m_resUpdates[LEFT][m_nbrEqs*iState+iVar] = m_divContFlxL[iState][iVar];
                m_resUpdates[RIGHT][m_nbrEqs*iState+iVar] = m_divContFlxR[iState][iVar];
              }
            }

            for (CFuint iSide = 0; iSide < 2; ++iSide)
            {
              // compute solution points Jacobian determinants
              m_solJacobDet[iSide] = m_cells[iSide]->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
            }

            // compute the diffusive face term contribution to the jacobian
            if ((*m_states[LEFT])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
            {
              computeBothJacobsDiffFaceTerm();
            }
            else if ((*m_states[LEFT])[0]->isParUpdatable())
            {
              computeOneJacobDiffFaceTerm(LEFT );
            }
            else if ((*m_states[RIGHT])[0]->isParUpdatable())
            {
              computeOneJacobDiffFaceTerm(RIGHT);
            }
          }
        }

        // release the cells
        m_cellBuilders[LEFT ]->releaseGE();
        m_cellBuilders[RIGHT]->releaseGE();

        m_cellFlags[cellIDL] = true;
        m_cellFlags[cellIDR] = true;
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionBlendingVS::computeBothJacobsDiffFaceTerm()
{
  // Get blending coefficient alpha
  DataHandle< CFreal > alphaData = socket_alpha.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices, proj vectors and make a backup of discontinuous fluxes
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[m_pertSide])[iSol]->getLocalID());
    }

    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_neighbCellFluxProjVects[m_pertSide][iDim] = m_cells[m_pertSide]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
    }

    // Loop over solution points to calculate the discontinuous flux.
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());

      prepareFluxComputation();

      // calculate the discontinuous flux projected on x, y, z-directions
      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxNeighb[m_pertSide][m_pertSol][iDim]);
        m_contFlxBackup[m_pertSide][m_pertSol][iDim] = m_contFlxNeighb[m_pertSide][m_pertSol][iDim];
      }
    }
  }

  // loop over left and right cell
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;

    // cell ID of the cell at the non-perturbed side
    const CFuint otherCellID = m_cells[iOtherSide]->getID();

    // term depending on iSide
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // term depending on iOtherSide
    const CFuint otherSideTerm = iOtherSide*m_nbrSolPnts;

    // Add the discontinuous gradient
    *m_cellStates = *(m_states[m_pertSide]);

    computeCellGradTerm(m_gradTermBefore);

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];

      // reset affected sol pnts
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_affectedSolPnts[LEFT][iSol] = false;
        m_affectedSolPnts[RIGHT][iSol] = false;
      }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

        // compute the perturbed gradients in the current cell
        computePerturbedGradientsAnalytical(m_pertSide);

        // compute the perturbed left and right states in the flx pnts
        computeFlxPntStatesAndGrads();

        // compute perturbed FI
        computeInterfaceFlxCorrection();

        // compute the perturbed corrections
        computeCorrection(m_pertSide, m_pertDivContFlx[m_pertSide]);
        computeCorrection(iOtherSide, m_pertDivContFlx[iOtherSide]);

        // put the perturbed and unperturbed corrections in the correct format
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_pertResUpdates[m_pertSide][m_nbrEqs*iState+iVar] = m_pertDivContFlx[m_pertSide][iState][iVar];
            m_pertResUpdates[iOtherSide][m_nbrEqs*iState+iVar] = m_pertDivContFlx[iOtherSide][iState][iVar];
          }
        }

        // compute the finite difference derivative of the face term (perturbed side rows)
        m_numJacob->computeDerivative(m_pertResUpdates[m_pertSide],m_resUpdates[m_pertSide],m_derivResUpdates);

        // multiply residual update derivatives with residual factor
        m_derivResUpdates *= resFactor;
        // SFV-6: no alpha scaling on diffusive Jacobian

        // add the derivative of the residual updates to the accumulator
        CFuint resUpdIdx = 0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
        {
          acc.addValues(iSol+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivResUpdates[resUpdIdx]);
        }

        // compute the perturbed diffusive residual in the other cell
        computePertCellDiffResiduals(iOtherSide);

        m_derivResUpdates = m_resUpdates[iOtherSide] + m_unpertAllCellDiffRes[otherCellID];

        // compute the finite difference derivative of the other cell the diffusive residual
        m_numJacob->computeDerivative(m_pertCellDiffRes,m_derivResUpdates,m_derivCellDiffRes);

        // multiply residual update derivatives with residual factor
        m_derivCellDiffRes *= resFactor;
        // SFV-6: no alpha scaling on diffusive Jacobian

        // add the derivative of the residual updates to the accumulator
        resUpdIdx = 0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
        {
          acc.addValues(iSol+otherSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivCellDiffRes[resUpdIdx]);
        }

        // Add internal cell contributions if needed
        const CFuint cellID = m_cells[m_pertSide]->getID();
        if (!m_cellFlags[cellID])
        {
          computeDivDiscontFlxNeighb(m_pertCellDiffRes,m_pertSide);

          // compute the finite difference derivative of the other cell the diffusive residual
          m_numJacob->computeDerivative(m_pertCellDiffRes,m_unpertAllCellDiffRes[cellID],m_derivCellDiffRes);

          // multiply residual update derivatives with residual factor
          m_derivCellDiffRes *= resFactor;
          // SFV-6: no alpha scaling on diffusive Jacobian

          // add the derivative of the residual updates to the accumulator
          resUpdIdx = 0;
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
          {
            acc.addValues(iSol+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivCellDiffRes[resUpdIdx]);
          }
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);

        // restore the gradients in the sol pnts
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            (*m_cellGrads[m_pertSide][iState])[iVar] = m_cellGradsBackUp[m_pertSide][iState][iVar];
            (*m_cellGrads[iOtherSide][iState])[iVar] = m_cellGradsBackUp[iOtherSide][iState][iVar];
          }

          for (CFuint iSide = 0; iSide < 2; ++iSide)
          {
            if (m_affectedSolPnts[iSide][iState])
            {
              // calculate the discontinuous flux projected on x, y, z-directions
              for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
              {
                m_contFlxNeighb[iSide][iState][iDim] = m_contFlxBackup[iSide][iState][iDim];
              }
            }
          }
        }
      }
    }
  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionBlendingVS::computeOneJacobDiffFaceTerm(const CFuint side)
{
  // Get blending coefficient alpha
  DataHandle< CFreal > alphaData = socket_alpha.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[m_pertSide])[iSol]->getLocalID());
    }

    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_neighbCellFluxProjVects[m_pertSide][iDim] = m_cells[m_pertSide]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
    }

    // Loop over solution points to calculate the discontinuous flux.
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());

      prepareFluxComputation();

      // calculate the discontinuous flux projected on x, y, z-directions
      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxNeighb[m_pertSide][m_pertSol][iDim]);
        m_contFlxBackup[m_pertSide][m_pertSol][iDim] = m_contFlxNeighb[m_pertSide][m_pertSol][iDim];
      }
    }
  }

  // loop over left and right cell
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;

    // cell ID of the cell at the non-perturbed side
    const CFuint otherCellID = m_cells[iOtherSide]->getID();

    // term depending on iSide
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // term depending on iOtherSide
    const CFuint otherSideTerm = iOtherSide*m_nbrSolPnts;

    // Add the discontinuous gradient
    *m_cellStates = *(m_states[m_pertSide]);

    computeCellGradTerm(m_gradTermBefore);

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];

      // reset affected sol pnts
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_affectedSolPnts[LEFT][iSol] = false;
        m_affectedSolPnts[RIGHT][iSol] = false;
      }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

        // compute the perturbed gradients in the current cell
        computePerturbedGradientsAnalytical(m_pertSide);

        // compute the perturbed left and right states in the flx pnts
        computeFlxPntStatesAndGrads();

        // compute perturbed FI
        computeInterfaceFlxCorrection();

        // compute the perturbed corrections
        computeCorrection(side, m_pertDivContFlx[side]);

        // put the perturbed and unperturbed corrections in the correct format
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_pertResUpdates[side][m_nbrEqs*iState+iVar] = m_pertDivContFlx[side][iState][iVar];
          }
        }

        // compute the finite difference derivative of the face term
        if (m_pertSide == side)
        {
          m_numJacob->computeDerivative(m_pertResUpdates[m_pertSide],m_resUpdates[m_pertSide],m_derivResUpdates);

          // multiply residual update derivatives with residual factor
          m_derivResUpdates *= resFactor;

          // SFV-6: no alpha scaling on diffusive Jacobian

          // add the derivative of the residual updates to the accumulator
          CFuint resUpdIdx = 0;
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
          {
            acc.addValues(iSol+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivResUpdates[resUpdIdx]);
          }

          // Add internal cell contributions if needed
          const CFuint cellID = m_cells[m_pertSide]->getID();
          if (!m_cellFlags[cellID])
          {
            computeDivDiscontFlxNeighb(m_pertCellDiffRes,m_pertSide);

            // compute the finite difference derivative of the other cell the diffusive residual
            m_numJacob->computeDerivative(m_pertCellDiffRes,m_unpertAllCellDiffRes[cellID],m_derivCellDiffRes);

            // multiply residual update derivatives with residual factor
            m_derivCellDiffRes *= resFactor;

            // SFV-6: no alpha scaling on diffusive Jacobian

            // add the derivative of the residual updates to the accumulator
            resUpdIdx = 0;
            for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
            {
              acc.addValues(iSol+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivCellDiffRes[resUpdIdx]);
            }
          }
        }
        else
        {
          // compute the perturbed diffusive residual in the other cell
          computePertCellDiffResiduals(iOtherSide);

          RealVector temp = m_resUpdates[iOtherSide] + m_unpertAllCellDiffRes[otherCellID];

          // compute the finite difference derivative of the other cell the diffusive residual
          m_numJacob->computeDerivative(m_pertCellDiffRes,temp,m_derivCellDiffRes);

          // multiply residual update derivatives with residual factor
          m_derivCellDiffRes *= resFactor;

          // SFV-6: no alpha scaling on diffusive Jacobian

          // add the derivative of the residual updates to the accumulator
          CFuint resUpdIdx = 0;
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
          {
            acc.addValues(iSol+otherSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivCellDiffRes[resUpdIdx]);
          }
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);

        // restore the gradients in the sol pnts
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            (*m_cellGrads[m_pertSide][iState])[iVar] = m_cellGradsBackUp[m_pertSide][iState][iVar];
            (*m_cellGrads[iOtherSide][iState])[iVar] = m_cellGradsBackUp[iOtherSide][iState][iVar];
          }

          for (CFuint iSide = 0; iSide < 2; ++iSide)
          {
            if (m_affectedSolPnts[iSide][iState])
            {
              // calculate the discontinuous flux projected on x, y, z-directions
              for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
              {
                m_contFlxNeighb[iSide][iState][iDim] = m_contFlxBackup[iSide][iState][iDim];
              }
            }
          }
        }
      }
    }
  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
