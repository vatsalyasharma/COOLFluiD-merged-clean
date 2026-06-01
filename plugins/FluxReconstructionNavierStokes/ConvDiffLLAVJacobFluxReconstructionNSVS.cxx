// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionNavierStokes/ConvDiffLLAVJacobFluxReconstructionNSVS.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvDiffLLAVJacobFluxReconstructionNSVS,
                       FluxReconstructionSolverData,
                       FluxReconstructionNavierStokesModule >
convDiffLLAVRHSJacobNSVSProvider("ConvDiffLLAVRHSJacobNSVS");

//////////////////////////////////////////////////////////////////////////////

ConvDiffLLAVJacobFluxReconstructionNSVS::ConvDiffLLAVJacobFluxReconstructionNSVS(const std::string& name) :
  ConvDiffLLAVJacobFluxReconstructionNSBaseVS(name),
  m_keepLLAVInJacob(false),
  m_diagValuesVS(),
  m_solPntsLocalCoordsVS(CFNULL)
{
  addConfigOptionsTo(this);

  m_keepLLAVInJacob = false;
  setParameter("KeepLLAVInJacob", &m_keepLLAVInJacob);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionNSVS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool, Config::DynamicOption<> >(
    "KeepLLAVInJacob",
    "When true AND Jacobian is frozen, replace stale matrix with LLAV+M/dt only. "
    "Set via .inter after epsilon and Jacobian are frozen.");
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionNSVS::configure(Config::ConfigArgs& args)
{
  ConvDiffLLAVJacobFluxReconstructionNSBaseVS::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionNSVS::setup()
{
  ConvDiffLLAVJacobFluxReconstructionNSBaseVS::setup();

  // allocate diagonal values vector
  m_diagValuesVS.resize(m_nbrSolPnts, 0.0);

  // get solution point local coordinates for Jacobian determinant
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  m_solPntsLocalCoordsVS = frLocalData[0]->getSolPntsLocalCoords();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionNSVS::execute()
{
  CFAUTOTRACE;

  // Step 1: call base class execute() — computes RHS always, Jacobian only when unfrozen
  ConvDiffLLAVJacobFluxReconstructionNSBaseVS::execute();

  // Step 2: check if we should do the LLAV-only second pass
  if (!m_keepLLAVInJacob) return;

  // Check freeze condition (same gate as base class, line 645)
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  const CFuint iterFreeze = getMethodData().getFreezeJacobIter();

  // If Jacobian is not frozen, the base class already computed the full Jacobian — nothing to do
  if (!getMethodData().freezeJacob()) return;
  if (iter < iterFreeze) return;
  const CFuint interval = iter - iterFreeze;
  if (interval % getMethodData().getFreezeJacobInterval() == 0) return;

  // We are in a frozen iteration and KeepLLAVInJacob is true.
  // Zero the PETSc matrix and recompute with LLAV + M/dt only.

  CFLog(INFO, "ConvDiffLLAVRHSJacobNSVS: rebuilding LLAV-only Jacobian + M/dt at iter " << iter << "\n");

  // Zero the global matrix
  SafePtr< LSSMatrix > jacobMatrix = m_lss->getMatrix();
  jacobMatrix->resetToZeroEntries();

  // Get cell volumes socket
  DataHandle< CFreal > cellVolumes = socket_cellVolumes.getDataHandle();

  // Get element type data
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // Get InnerCells and InnerFaces TRS
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // Get face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size() - 1;

  // Set up geo builders
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;

  CellToFaceGEBuilder::GeoData& geoDataCBL = m_cellBuilders[LEFT]->getDataGE();
  geoDataCBL.trs = cells;
  CellToFaceGEBuilder::GeoData& geoDataCBR = m_cellBuilders[RIGHT]->getDataGE();
  geoDataCBR.trs = cells;

  // Reset cell flags for the second pass
  for (CFuint iCell = 0; iCell < m_cellFlags.size(); ++iCell)
  {
    m_cellFlags[iCell] = false;
  }

  // Get sol pnt normals and volumes
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
  DataHandle< CFreal > volumes = socket_volumes.getDataHandle();

  //// Second face loop: compute LLAV Jacobian + M/dt
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient + 1];

    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      m_faceID = faceID;

      // Build face
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();

      // Get neighboring cells
      m_cells[LEFT]  = m_face->getNeighborGeo(LEFT);
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      m_states[LEFT]  = m_cells[LEFT]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();

      // Cell volumes
      m_cellVolume[LEFT]  = cellVolumes[m_cells[LEFT]->getID()];
      m_cellVolume[RIGHT] = cellVolumes[m_cells[RIGHT]->getID()];

      // Build cell GEs
      const CFuint cellIDL = m_face->getNeighborGeo(LEFT)->getID();
      geoDataCBL.idx = cellIDL;
      m_cells[LEFT] = m_cellBuilders[LEFT]->buildGE();

      const CFuint cellIDR = m_face->getNeighborGeo(RIGHT)->getID();
      geoDataCBR.idx = cellIDR;
      m_cells[RIGHT] = m_cellBuilders[RIGHT]->buildGE();

      // Set face data (populates m_epsilonLR, m_cellGrads, m_cellGradsAV, normals, etc.)
      setFaceData(m_face->getID());

      // Compute sol pnt Jacobian determinants and sol pnt epsilons
      for (CFuint iSide = 0; iSide < 2; ++iSide)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          m_solJacobDet[iSide][iSol] = volumes[(*(m_states[iSide]))[iSol]->getLocalID()];

          // Interpolate epsilon from nodes to sol pnts
          m_solEpsilons[iSide][iSol] = 0.0;
          m_cellNodes[iSide] = m_cells[iSide]->getNodes();
          for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
          {
            const CFuint nodeIdx = (*(m_cellNodes[iSide]))[iNode]->getLocalID();
            m_solEpsilons[iSide][iSol] += m_nodePolyValsAtSolPnts[iSol][iNode] *
                                           m_nodeEpsilons[nodeIdx] / m_nbNodeNeighbors[nodeIdx];
          }

          // Sol pnt normals (flux projection vectors)
          const CFuint solID = (*(m_states[iSide]))[iSol]->getLocalID();
          for (CFuint iDim = 0; iDim < m_dim + m_ndimplus; ++iDim)
          {
            for (CFuint jDim = 0; jDim < m_dim; ++jDim)
            {
              m_neighbCellFluxProjVects[iSide][iDim][iSol][jDim] =
                solPntNormals[solID * (m_dim + m_ndimplus) * m_dim + iDim * m_dim + jDim];
            }
          }
        }
      }

      // Compute LLAV Jacobian for this face pair
      if ((*m_states[LEFT])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
      {
        computeLLAVOnlyBothJacobs();
      }
      else if ((*m_states[LEFT])[0]->isParUpdatable())
      {
        computeLLAVOnlyOneJacob(LEFT);
      }
      else if ((*m_states[RIGHT])[0]->isParUpdatable())
      {
        computeLLAVOnlyOneJacob(RIGHT);
      }

      // Add M/dt diagonal for cells not yet processed
      for (CFuint iSide = 0; iSide < 2; ++iSide)
      {
        const CFuint cellID = m_cells[iSide]->getID();
        if (!m_cellFlags[cellID] && (*(m_states[iSide]))[0]->isParUpdatable())
        {
          addTimeDiagonal(iSide);
        }
      }

      // Release cells
      m_cellBuilders[LEFT]->releaseGE();
      m_cellBuilders[RIGHT]->releaseGE();

      m_cellFlags[cellIDL] = true;
      m_cellFlags[cellIDR] = true;

      // Release face
      m_faceBuilder->releaseGE();
    }
  }

  // Finalize matrix assembly
  jacobMatrix->finalAssembly();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionNSVS::computeLLAVOnlyBothJacobs()
{
  CFLog(VERBOSE, "computeLLAVOnlyBothJacobs\n");

  // Dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // Set row/col indices for both sides
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx, (*m_states[m_pertSide])[iSol]->getLocalID());
    }
  }

  // Compute analytical gradient-to-state Jacobian (needed for LLAV terms)
  computeGradToStateJacobianAna();

  //// Add LLAV Jacobian contributions (LLAV-only: skip m_fluxJacobian, m_riemannFluxJacobian, etc.)

  // Cell (discontinuous) LLAV part
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    if (!m_cellFlags[m_cells[m_pertSide]->getID()])
    {
      const CFuint pertSideTerm = m_pertSide * m_nbrSolPnts;

      for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
      {
        // Compute dCons/dState Jacobian at this perturbation point via FD
        {
          State* statePtr = (*m_states[m_pertSide])[m_pertSol];
          m_origConsState = *m_updateToSolutionVecTrans->transform(statePtr);
          const CFreal eps_fd = 1.0e-8;
          for (CFuint j = 0; j < m_nbrEqs; ++j)
          {
            const CFreal origVal = (*statePtr)[j];
            (*statePtr)[j] += eps_fd;
            const State* pertConsPtr = m_updateToSolutionVecTrans->transform(statePtr);
            for (CFuint k = 0; k < m_nbrEqs; ++k)
            {
              m_avConsToStateJac(k, j) = ((*pertConsPtr)[k] - m_origConsState[k]) / eps_fd;
            }
            (*statePtr)[j] = origVal;
          }
        }

        for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
        {
          // --- Cell volume LLAV (discontinuous gradient part) ---
          // First part: derivative through sol pnt polynomial derivatives
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
          {
            const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];

            for (CFuint jSol = 0; jSol < m_nbrSolSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_solSolDep)[kSolIdx][jSol];

              m_tempFlux = 0.0;

              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal dl = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdx];

                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal dl_dqdu = dl * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];

                  CFreal llavPart = m_solEpsilons[m_pertSide][kSolIdx] *
                                    m_neighbCellFluxProjVects[m_pertSide][iDim][kSolIdx][jDim];
                  llavPart *= dl_dqdu;
                  addLLAVJacobToFlux(llavPart);
                }
              }

              acc.addValues(jSolIdx + pertSideTerm, m_pertSol + pertSideTerm, m_pertVar, &m_tempFlux[0]);
            }
          }

          // Second part: through flux-point correction function
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
          {
            const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];

            for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
            {
              const CFuint flxIdx = (*m_solFlxDep)[kSolIdx][iFlxPnt];
              const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
              const CFreal l = (*m_solPolyValsAtFlxPnts)[flxIdx][kSolIdx];

              m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
              for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
              {
                const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];

                m_tempFlux = 0.0;

                const CFreal divh_l = -m_corrFctDiv[jSolIdx][flxIdx] * l;

                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqdu = divh_l *
                    m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];

                  CFreal llavPart = divh_l_dqdu * m_solEpsilons[m_pertSide][kSolIdx];
                  llavPart *= m_neighbCellFluxProjVects[m_pertSide][dim][kSolIdx][jDim];
                  addLLAVJacobToFlux(llavPart);
                }

                acc.addValues(jSolIdx + pertSideTerm, m_pertSol + pertSideTerm, m_pertVar, &m_tempFlux[0]);
              }
            }
          }

        } // m_pertVar
      } // m_pertSol
    } // cellFlags
  } // m_pertSide (cell part)

  // Face (Riemann) LLAV part
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    const CFuint pertSideTerm = m_pertSide * m_nbrSolPnts;
    const CFuint otherSideTerm = iOtherSide * m_nbrSolPnts;

    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint flxPntIdxThis  = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
        const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlxPnt];

        const CFreal halfFaceJacob = 0.5 * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        const CFreal epsilon = 0.5 * (m_epsilonLR[LEFT][iFlxPnt] + m_epsilonLR[RIGHT][iFlxPnt]);

        m_nbrSolDep = ((*m_flxSolDep)[flxPntIdxThis]).size();

        // --- Same-cell Riemann gradient LLAV ---
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
        {
          const CFuint jSolIdxThis  = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

          CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
          const CFreal divh_halfFaceJacob = divh * halfFaceJacob;

          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];

            if (m_addRiemannToGradJacob)
            {
              // This side
              m_tempFlux = 0.0;
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx      = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
                const CFuint kSolIdxOther  = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];

                const CFreal divh_halfFaceJacob_l      = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
                const CFreal divh_halfFaceJacob_lOther = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];

                for (CFuint iDim = 0; iDim < m_dim; ++iDim)
                {
                  const CFreal divh_halfFaceJacob_l_dqdu =
                    divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                  const CFreal divh_halfFaceJacob_l_dqduOther =
                    divh_halfFaceJacob_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];

                  const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                  addLLAVJacobToFlux(llavPart * divh_halfFaceJacob_l_dqdu);
                  addLLAVJacobToFlux(llavPart * divh_halfFaceJacob_l_dqduOther);
                }
              }
              acc.addValues(jSolIdxThis + pertSideTerm, pertSolIdx + pertSideTerm, m_pertVar, &m_tempFlux[0]);

              // Other side (cross-cell)
              divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
              const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

              m_tempFlux = 0.0;
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx      = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
                const CFuint kSolIdxOther  = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];

                const CFreal divh_halfFaceJacobOther_lThis  = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
                const CFreal divh_halfFaceJacobOther_lOther = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];

                for (CFuint iDim = 0; iDim < m_dim; ++iDim)
                {
                  const CFreal divh_halfFaceJacobOther_lThis_dqduThis =
                    divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                  const CFreal divh_halfFaceJacobOther_lOther_dqduOther =
                    divh_halfFaceJacobOther_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];

                  const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                  addLLAVJacobToFlux(llavPart * divh_halfFaceJacobOther_lThis_dqduThis);
                  addLLAVJacobToFlux(llavPart * divh_halfFaceJacobOther_lOther_dqduOther);
                }
              }
              acc.addValues(jSolIdxOther + otherSideTerm, pertSolIdx + pertSideTerm, m_pertVar, &m_tempFlux[0]);
            } // addRiemannToGradJacob

          } // m_pertSol (same-cell Riemann)

          // Cross-cell Riemann gradient LLAV (for sol pnts not on the flux stencil)
          if (m_addRiemannToGradCrossCellJacob)
          {
            // Mark which sol pnts are already covered
            for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
            {
              m_needToAddSolPnt[iSol] = true;
            }
            for (CFuint iPS = 0; iPS < m_nbrSolDep; ++iPS)
            {
              m_needToAddSolPnt[(*m_flxSolDep)[flxPntIdxThis][iPS]] = false;
            }

            for (CFuint pertSolIdx = 0; pertSolIdx < m_nbrSolPnts; ++pertSolIdx)
            {
              if (!m_needToAddSolPnt[pertSolIdx]) continue;

              CFuint dependingKSol = 1000;
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
                for (CFuint lSol = 0; lSol < m_nbrSolSolDep; ++lSol)
                {
                  if ((*m_solSolDep)[pertSolIdx][lSol] == kSolIdx)
                  {
                    dependingKSol = kSolIdx;
                    break;
                  }
                }
              }
              if (dependingKSol == 1000) continue;

              const CFreal divh_hFJ_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][dependingKSol];

              // This side
              m_tempFlux = 0.0;
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal dqdu = divh_hFJ_l *
                  m_gradientStateJacobian[m_pertSide][dependingKSol][m_pertSide][pertSolIdx][iDim];
                const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                addLLAVJacobToFlux(llavPart * dqdu);
              }
              acc.addValues(jSolIdxThis + pertSideTerm, pertSolIdx + pertSideTerm, m_pertVar, &m_tempFlux[0]);

              // Other side
              divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
              const CFreal divh_hFJO_lThis = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide]
                                           * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][dependingKSol];
              m_tempFlux = 0.0;
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal dqdu = divh_hFJO_lThis *
                  m_gradientStateJacobian[m_pertSide][dependingKSol][m_pertSide][pertSolIdx][iDim];
                const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                addLLAVJacobToFlux(llavPart * dqdu);
              }
              acc.addValues(jSolIdxOther + otherSideTerm, pertSolIdx + pertSideTerm, m_pertVar, &m_tempFlux[0]);
            } // pertSolIdx
          } // addRiemannToGradCrossCellJacob

        } // jSolPnt

        //// Cross-element gradient part (flux derivative through other-side gradient)
        if (m_addFluxToGradCrossCellJacob)
        {
          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];

            // Through correction function on other side
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];

              for (CFuint iInfluencedFlx = 0; iInfluencedFlx < m_nbrFlxDep; ++iInfluencedFlx)
              {
                const CFuint iInfluencedFlxIdx = (*m_solFlxDep)[kSolIdxOther][iInfluencedFlx];
                const CFuint dimOther = (*m_flxPntFlxDim)[iInfluencedFlxIdx];
                const CFreal lOther = (*m_solPolyValsAtFlxPnts)[iInfluencedFlxIdx][kSolIdxOther];

                for (CFuint jSolPnt2 = 0; jSolPnt2 < m_nbrSolSolDep; ++jSolPnt2)
                {
                  const CFuint jSolIdx2 = (*m_solSolDep)[kSolIdxOther][jSolPnt2];

                  m_tempFlux = 0.0;
                  const CFreal divh_lOther = -m_corrFctDiv[jSolIdx2][iInfluencedFlxIdx] * lOther;

                  for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                  {
                    const CFreal divh_l_dqduOther = divh_lOther *
                      m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];

                    CFreal llavPart = divh_l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
                    llavPart *= m_neighbCellFluxProjVects[iOtherSide][dimOther][kSolIdxOther][jDim];
                    addLLAVJacobToFlux(llavPart);
                  }

                  acc.addValues(jSolIdx2 + otherSideTerm, pertSolIdx + pertSideTerm, m_pertVar, &m_tempFlux[0]);
                }
              }
            }

            // Through sol pnt derivative on other side
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];

              for (CFuint jSolPnt2 = 0; jSolPnt2 < m_nbrSolSolDep; ++jSolPnt2)
              {
                const CFuint jSolIdx2 = (*m_solSolDep)[kSolIdxOther][jSolPnt2];

                m_tempFlux = 0.0;

                for (CFuint iDim = 0; iDim < m_dim; ++iDim)
                {
                  const CFreal lOther = (*m_solPolyDerivAtSolPnts)[jSolIdx2][iDim][kSolIdxOther];

                  for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                  {
                    const CFreal l_dqduOther = lOther *
                      m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];

                    CFreal llavPart = l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
                    llavPart *= m_neighbCellFluxProjVects[iOtherSide][iDim][kSolIdxOther][jDim];
                    addLLAVJacobToFlux(llavPart);
                  }
                }

                acc.addValues(jSolIdx2 + otherSideTerm, pertSolIdx + pertSideTerm, m_pertVar, &m_tempFlux[0]);
              }
            }

          } // m_pertSol (cross-cell flux gradient)
        } // addFluxToGradCrossCellJacob

      } // iFlxPnt
    } // m_pertVar
  } // m_pertSide (face part)

  if (getMethodData().doComputeJacobian())
  {
    m_lss->getMatrix()->addValues(acc);
  }

  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionNSVS::computeLLAVOnlyOneJacob(const CFuint side)
{
  CFLog(VERBOSE, "computeLLAVOnlyOneJacob, side = " << side << "\n");

  // This is structurally identical to computeLLAVOnlyBothJacobs.
  // The only difference is that on MPI partition boundaries, only the
  // "side" cell is par-updatable. The full two-sided stencil is still
  // needed for cross-cell gradient coupling, so we assemble the same
  // terms but only add to the matrix for the updatable side.

  BlockAccumulator& acc = *m_acc;

  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx, (*m_states[m_pertSide])[iSol]->getLocalID());
    }
  }

  computeGradToStateJacobianAna();

  // Cell LLAV part — only for the updatable side
  m_pertSide = side;
  if (!m_cellFlags[m_cells[m_pertSide]->getID()])
  {
    const CFuint pertSideTerm = m_pertSide * m_nbrSolPnts;

    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dCons/dState Jacobian
      {
        State* statePtr = (*m_states[m_pertSide])[m_pertSol];
        m_origConsState = *m_updateToSolutionVecTrans->transform(statePtr);
        const CFreal eps_fd = 1.0e-8;
        for (CFuint j = 0; j < m_nbrEqs; ++j)
        {
          const CFreal origVal = (*statePtr)[j];
          (*statePtr)[j] += eps_fd;
          const State* pertConsPtr = m_updateToSolutionVecTrans->transform(statePtr);
          for (CFuint k = 0; k < m_nbrEqs; ++k)
          {
            m_avConsToStateJac(k, j) = ((*pertConsPtr)[k] - m_origConsState[k]) / eps_fd;
          }
          (*statePtr)[j] = origVal;
        }
      }

      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // First part: through sol pnt poly derivatives
        for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
        {
          const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];

          for (CFuint jSol = 0; jSol < m_nbrSolSolDep; ++jSol)
          {
            const CFuint jSolIdx = (*m_solSolDep)[kSolIdx][jSol];

            m_tempFlux = 0.0;

            for (CFuint iDim = 0; iDim < m_dim; ++iDim)
            {
              const CFreal dl = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdx];

              for (CFuint jDim = 0; jDim < m_dim; ++jDim)
              {
                const CFreal dl_dqdu = dl * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                CFreal llavPart = m_solEpsilons[m_pertSide][kSolIdx] *
                                  m_neighbCellFluxProjVects[m_pertSide][iDim][kSolIdx][jDim];
                llavPart *= dl_dqdu;
                addLLAVJacobToFlux(llavPart);
              }
            }

            acc.addValues(jSolIdx + pertSideTerm, m_pertSol + pertSideTerm, m_pertVar, &m_tempFlux[0]);
          }
        }

        // Second part: through correction function
        for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
        {
          const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];

          for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
          {
            const CFuint flxIdx = (*m_solFlxDep)[kSolIdx][iFlxPnt];
            const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
            const CFreal l = (*m_solPolyValsAtFlxPnts)[flxIdx][kSolIdx];

            m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];

              m_tempFlux = 0.0;
              const CFreal divh_l = -m_corrFctDiv[jSolIdx][flxIdx] * l;

              for (CFuint jDim = 0; jDim < m_dim; ++jDim)
              {
                const CFreal divh_l_dqdu = divh_l *
                  m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                CFreal llavPart = divh_l_dqdu * m_solEpsilons[m_pertSide][kSolIdx];
                llavPart *= m_neighbCellFluxProjVects[m_pertSide][dim][kSolIdx][jDim];
                addLLAVJacobToFlux(llavPart);
              }

              acc.addValues(jSolIdx + pertSideTerm, m_pertSol + pertSideTerm, m_pertVar, &m_tempFlux[0]);
            }
          }
        }
      } // m_pertVar
    } // m_pertSol
  } // cellFlags

  // Face LLAV part — same as BothJacobs (both pertSides needed for cross-cell coupling)
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    const CFuint pertSideTerm = m_pertSide * m_nbrSolPnts;
    const CFuint otherSideTerm = iOtherSide * m_nbrSolPnts;

    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint flxPntIdxThis  = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
        const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlxPnt];

        const CFreal halfFaceJacob = 0.5 * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        const CFreal epsilon = 0.5 * (m_epsilonLR[LEFT][iFlxPnt] + m_epsilonLR[RIGHT][iFlxPnt]);

        m_nbrSolDep = ((*m_flxSolDep)[flxPntIdxThis]).size();

        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
        {
          const CFuint jSolIdxThis  = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

          CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
          const CFreal divh_halfFaceJacob = divh * halfFaceJacob;

          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];

            if (m_addRiemannToGradJacob)
            {
              // This side
              m_tempFlux = 0.0;
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx      = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
                const CFuint kSolIdxOther  = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];

                const CFreal dhfl      = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
                const CFreal dhflOther = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];

                for (CFuint iDim = 0; iDim < m_dim; ++iDim)
                {
                  const CFreal dhfl_dqdu      = dhfl * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                  const CFreal dhfl_dqduOther = dhflOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];

                  const CFreal lp = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                  addLLAVJacobToFlux(lp * dhfl_dqdu);
                  addLLAVJacobToFlux(lp * dhfl_dqduOther);
                }
              }
              acc.addValues(jSolIdxThis + pertSideTerm, pertSolIdx + pertSideTerm, m_pertVar, &m_tempFlux[0]);

              // Other side
              divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
              const CFreal dhfjo = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

              m_tempFlux = 0.0;
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx      = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
                const CFuint kSolIdxOther  = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];

                const CFreal dhfjo_lThis  = dhfjo * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
                const CFreal dhfjo_lOther = dhfjo * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];

                for (CFuint iDim = 0; iDim < m_dim; ++iDim)
                {
                  const CFreal dqduThis  = dhfjo_lThis * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                  const CFreal dqduOther = dhfjo_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];

                  const CFreal lp = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                  addLLAVJacobToFlux(lp * dqduThis);
                  addLLAVJacobToFlux(lp * dqduOther);
                }
              }
              acc.addValues(jSolIdxOther + otherSideTerm, pertSolIdx + pertSideTerm, m_pertVar, &m_tempFlux[0]);
            }
          } // m_pertSol
        } // jSolPnt
      } // iFlxPnt
    } // m_pertVar
  } // m_pertSide

  if (getMethodData().doComputeJacobian())
  {
    m_lss->getMatrix()->addValues(acc);
  }

  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionNSVS::addTimeDiagonal(const CFuint side)
{
  // This adds M/dt diagonal to the Jacobian for the given cell side.
  // Pattern from PseudoSteadyStdTimeRHSJacob (pseudo-steady case).

  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  const CFreal resFactor = getMethodData().getResFactor();

  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();

  const CFreal invCellVolume = 1.0 / m_cellVolume[side];
  const CFuint firstStateID = (*(m_states[side]))[0]->getLocalID();
  const CFreal updateCoeffDivCFL = invCellVolume * updateCoeff[firstStateID] / cfl;

  // Use single-cell accumulator
  BlockAccumulator& accSC = *m_accSC;

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    accSC.setRowColIndex(iSol, (*(m_states[side]))[iSol]->getLocalID());
  }

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFreal diagValue = updateCoeffDivCFL * m_solJacobDet[side][iSol] * resFactor;

    State* currState = (*(m_states[side]))[iSol];

    // Compute dSol/dUpdate via FD (same pattern as PseudoSteadyStdTimeRHSJacob)
    m_tempSolVarState = *m_updateToSolutionVecTrans->transform(currState);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // Perturb
      const CFreal origVal = (*currState)[iEq];
      const CFreal eps_fd = 1.0e-8;
      (*currState)[iEq] += eps_fd;

      const RealVector& pertState =
        static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(currState));

      // dSol/dUpdate column iEq
      for (CFuint kEq = 0; kEq < m_nbrEqs; ++kEq)
      {
        m_tempFlux[kEq] = (pertState[kEq] - m_tempSolVarState[kEq]) / eps_fd * diagValue;
      }

      accSC.addValues(iSol, iSol, iEq, &m_tempFlux[0]);

      // Restore
      (*currState)[iEq] = origVal;
    }
  }

  if (getMethodData().doComputeJacobian())
  {
    m_lss->getMatrix()->addValues(accSC);
  }

  accSC.reset();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
