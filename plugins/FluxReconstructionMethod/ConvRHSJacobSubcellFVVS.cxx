// Copyright (C) 2026
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BlockAccumulator.hh"
#include "Framework/CFSide.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "FluxReconstructionMethod/ConvRHSJacobSubcellFVVS.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
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

MethodCommandProvider< ConvRHSJacobSubcellFVVS,
                       FluxReconstructionSolverData,
                       FluxReconstructionModule >
  ConvRHSJacobSubcellFVVSProvider("ConvRHSJacobSubcellFVVS");

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("AlphaEps",
    "Skip subcell FV when alpha <= AlphaEps (default 1e-12).");
}

//////////////////////////////////////////////////////////////////////////////

ConvRHSJacobSubcellFVVS::ConvRHSJacobSubcellFVVS(const std::string& name) :
  ConvRHSJacobFluxReconstruction(name),
  socket_alpha("alpha"),
  m_solPnts1D(),
  m_subcellInterfaces1D(),
  m_invSubcellWidth1D(),
  m_nbrSolPnts1D(0),
  m_tensorBased(false),
  m_subcellIfaceCoords(),
  m_subcellIfaceDimList(),
  m_subcellIfaceLeftIdx(),
  m_subcellIfaceRightIdx(),
  m_subcellIfaceInvWidthL(),
  m_subcellIfaceInvWidthR(),
  m_subcellFVRes(),
  m_alphaEps(1.0e-12)
{
  addConfigOptionsTo(this);
  setParameter("AlphaEps", &m_alphaEps);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::configure ( Config::ConfigArgs& args )
{
  ConvRHSJacobFluxReconstruction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
ConvRHSJacobSubcellFVVS::needsSockets()
{
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > result =
    ConvRHSJacobFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::setup()
{
  CFAUTOTRACE;

  ConvRHSJacobFluxReconstruction::setup();

  // resize subcell FV residual storage
  m_subcellFVRes.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_subcellFVRes[iSol].resize(m_nbrEqs);
  }

  // unit normal for subcell FV Riemann flux
  m_subcellUnitNormal.resize(m_dim);

  // --- AJ-1: Analytical subcell FV Jacobian storage ---
  // Determine equation structure: species first, then Euler (momentum + energy)
  // For CNEQ 2D: [ρ_N, ρ_N2, ρu, ρv, ρE], nbSpecies=2, dim=2
  m_nbSpecies = m_nbrEqs - m_dim - 1;
  m_eulerStartIdx = m_nbSpecies;

  // Physical data size and first scalar variable index
  // EulerTerm: {RHO=0,..,GAMMA=10,XP=11,YP=12,ZP=13}, firstScalarVar=14
  m_firstScalarPdata = 14;

  m_pdataSize = m_pData.size();

  // Per-solution-point storage
  m_pdataPerSol.resize(m_nbrSolPnts);
  m_dPdUPerSol.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_pdataPerSol[iSol].resize(m_pdataSize, 0.0);
    m_dPdUPerSol[iSol].resize(m_pdataSize, m_nbrEqs, 0.0);
  }
  m_pdataPerturbed.resize(m_pdataSize, 0.0);

  // Flux Jacobian blocks
  m_dFluxdUL.resize(m_nbrEqs, m_nbrEqs, 0.0);
  m_dFluxdUR.resize(m_nbrEqs, m_nbrEqs, 0.0);

  // default: assume not tensor-based
  m_tensorBased = false;

  // get 1D solution point coordinates
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  if (frLocalData.empty())
  {
    return;
  }

  SafePtr< vector< CFreal > > sol1D = frLocalData[0]->getSolPntsLocalCoord1D();
  if (sol1D.isNull() || sol1D->empty())
  {
    return;
  }

  m_solPnts1D = *sol1D;
  m_nbrSolPnts1D = static_cast<CFuint>(m_solPnts1D.size());

  // Check tensor-product compatibility: N1D^m_dim must equal nbrSolPnts
  // Supports QUAD (2D) and HEXA (3D). Non-tensor-product elements
  // (TRIAG, TETRA, PRISM) will fail this check and fall back to pure FR.
  CFuint tensorPnts = 1;
  for (CFuint d = 0; d < m_dim; ++d) tensorPnts *= m_nbrSolPnts1D;

  if (tensorPnts != m_nbrSolPnts)
  {
    CFLog(WARN, "ConvRHSJacobSubcellFVVS: non tensor-product element ("
          << m_nbrSolPnts << " sol pnts, expected " << tensorPnts
          << " for " << m_dim << "D). Subcell FV disabled.\n");
    return;
  }

  m_tensorBased = true;

  // build 1D subcell interfaces (shared by all directions)
  m_subcellInterfaces1D.resize(m_nbrSolPnts1D + 1);
  m_subcellInterfaces1D[0] = -1.0;
  for (CFuint i = 0; i < m_nbrSolPnts1D - 1; ++i)
  {
    m_subcellInterfaces1D[i + 1] = 0.5 * (m_solPnts1D[i] + m_solPnts1D[i + 1]);
  }
  m_subcellInterfaces1D[m_nbrSolPnts1D] = 1.0;

  m_invSubcellWidth1D.resize(m_nbrSolPnts1D);
  for (CFuint i = 0; i < m_nbrSolPnts1D; ++i)
  {
    const CFreal width = m_subcellInterfaces1D[i + 1] - m_subcellInterfaces1D[i];
    m_invSubcellWidth1D[i] = (width > 0.0) ? (1.0 / width) : 0.0;
  }

  // Compute strides for tensor-product indexing:
  // idx(i_0, i_1, ..., i_{D-1}) = sum_d i_d * stride[d]
  // where i_0 = ksi (slowest), i_{D-1} = last dir (fastest)
  // stride[D-1] = 1, stride[d] = stride[d+1] * N1D
  std::vector< CFuint > strides(m_dim);
  strides[m_dim - 1] = 1;
  for (CFint d = static_cast<CFint>(m_dim) - 2; d >= 0; --d)
  {
    strides[d] = strides[d + 1] * m_nbrSolPnts1D;
  }

  // Number of solution points on each subcell interface face = N1D^(D-1)
  const CFuint nPtsOnFace = m_nbrSolPnts / m_nbrSolPnts1D;

  // Build subcell interface lists for each direction (generic D-dimensional)
  m_subcellIfaceCoords.resize(m_dim);
  m_subcellIfaceDimList.resize(m_dim);
  m_subcellIfaceLeftIdx.resize(m_dim);
  m_subcellIfaceRightIdx.resize(m_dim);
  m_subcellIfaceInvWidthL.resize(m_dim);
  m_subcellIfaceInvWidthR.resize(m_dim);

  for (CFuint dir = 0; dir < m_dim; ++dir)
  {
    m_subcellIfaceCoords[dir].clear();
    m_subcellIfaceDimList[dir].clear();
    m_subcellIfaceLeftIdx[dir].clear();
    m_subcellIfaceRightIdx[dir].clear();
    m_subcellIfaceInvWidthL[dir].clear();
    m_subcellIfaceInvWidthR[dir].clear();

    // Collect non-dir dimensions and their decomposition strides
    std::vector< CFuint > otherDims;
    std::vector< CFuint > otherStrides;
    for (CFuint d = 0; d < m_dim; ++d)
    {
      if (d != dir) otherDims.push_back(d);
    }
    if (!otherDims.empty())
    {
      otherStrides.resize(otherDims.size());
      otherStrides.back() = 1;
      for (CFint od = static_cast<CFint>(otherDims.size()) - 2; od >= 0; --od)
      {
        otherStrides[od] = otherStrides[od + 1] * m_nbrSolPnts1D;
      }
    }

    // Loop over interior interfaces in direction dir (N1D - 1 interfaces)
    for (CFuint iface = 0; iface < m_nbrSolPnts1D - 1; ++iface)
    {
      const CFreal mid = 0.5 * (m_solPnts1D[iface] + m_solPnts1D[iface + 1]);

      // Loop over all points on this interface plane
      for (CFuint ip = 0; ip < nPtsOnFace; ++ip)
      {
        RealVector coord(m_dim);
        coord[dir] = mid;

        // Decompose ip into indices for non-dir dimensions
        CFuint idxBase = iface * strides[dir];
        CFuint rem = ip;
        for (CFuint od = 0; od < otherDims.size(); ++od)
        {
          const CFuint d = otherDims[od];
          const CFuint idx_d = rem / otherStrides[od];
          rem %= otherStrides[od];
          coord[d] = m_solPnts1D[idx_d];
          idxBase += idx_d * strides[d];
        }

        m_subcellIfaceCoords[dir].push_back(coord);
        m_subcellIfaceDimList[dir].push_back(dir);
        m_subcellIfaceLeftIdx[dir].push_back(idxBase);
        m_subcellIfaceRightIdx[dir].push_back(idxBase + strides[dir]);
        m_subcellIfaceInvWidthL[dir].push_back(m_invSubcellWidth1D[iface]);
        m_subcellIfaceInvWidthR[dir].push_back(m_invSubcellWidth1D[iface + 1]);
      }
    }
  }

  // SFV-5: Build boundary subcell mapping for each flux point.
  // Maps each face flux point to the adjacent boundary subcell (solution point)
  // and its inverse subcell width in the face-normal direction.
  const CFuint nbrFlxPnts = m_flxPntsLocalCoords->size();
  m_flxPntBndSolIdx.resize(nbrFlxPnts);
  m_flxPntBndInvWidth.resize(nbrFlxPnts);

  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    const RealVector& flxCoord = (*m_flxPntsLocalCoords)[iFlx];
    const CFuint flxDim = (*m_flxPntFlxDim)[iFlx];

    // Boundary subcell index: at +1 boundary → last subcell, at -1 → first
    const CFuint bndIdx1D = (flxCoord[flxDim] > 0.0) ? (m_nbrSolPnts1D - 1) : 0;

    // Compute tensor-product solution point index
    CFuint solIdx = bndIdx1D * strides[flxDim];
    for (CFuint d = 0; d < m_dim; ++d)
    {
      if (d != flxDim)
      {
        // Find nearest solution point in transverse direction
        CFuint nearestIdx = 0;
        CFreal minDist = std::abs(flxCoord[d] - m_solPnts1D[0]);
        for (CFuint i = 1; i < m_nbrSolPnts1D; ++i)
        {
          const CFreal dist = std::abs(flxCoord[d] - m_solPnts1D[i]);
          if (dist < minDist)
          {
            minDist = dist;
            nearestIdx = i;
          }
        }
        solIdx += nearestIdx * strides[d];
      }
    }

    m_flxPntBndSolIdx[iFlx] = solIdx;
    m_flxPntBndInvWidth[iFlx] = m_invSubcellWidth1D[bndIdx1D];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::unsetup()
{
  ConvRHSJacobFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "ConvRHSJacobSubcellFVVS::execute()\n");

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;

  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size() - 1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;

  //// Loop over faces to calculate fluxes and interface fluxes in the flux points

  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    CFLog(VERBOSE, "Orient = " << m_orient << "\n");
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient + 1];

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

      // if one of the neighbouring cells is parallel updatable or if the gradients have to be computed, set the bnd face data
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable() || hasDiffTerm)
      {
        // set the bnd face data
        setFaceData(m_face->getID());

        // compute the states in the flx pnts
        computeFlxPntStates();

        // compute the interface flux
        computeInterfaceFlxCorrection();

        // compute the wave speed updates
        computeWaveSpeedUpdates(m_waveSpeedUpd);

        // update the wave speed
        updateWaveSpeed();
      }

      // if one of the neighbouring cells is parallel updatable, compute the correction flux
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
        // compute the correction for the left neighbour
        computeCorrection(LEFT, m_divContFlxL);

        // compute the correction for the right neighbour
        computeCorrection(RIGHT, m_divContFlxR);

        // NOTE: Face corrections are NOT blended here. The divhFD term in
        // computeDivDiscontFlx IS scaled by (1-alpha) (SFV-5), which removes
        // the extrapolated flux contribution at high alpha, leaving only the
        // Riemann flux (via the face correction) to provide boundary coupling.
        // Full face loop blending (scaling corrections by (1-alpha) + FV
        // boundary contributions) was tested but causes instability because
        // invSubcellWidth and corrFctDiv have different magnitudes, creating
        // a large perturbation when alpha transitions from 0 to nonzero.

        // update RHS
        updateRHSBothSides();
      }

      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
        // compute the face correction term of the corrected gradients
        computeGradientFaceCorrections();
      }

      // compute the contribution to the face numerical jacobian (Fix FJ1: gate with freeze check)
      {
        const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
        const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
        const CFuint interval = iter - iterFreeze;
        if (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0)
        {
          if ((*m_states[LEFT])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
          {
            computeBothJacobs();
          }
          else if ((*m_states[LEFT])[0]->isParUpdatable())
          {
            computeOneJacob(LEFT);
          }
          else if ((*m_states[RIGHT])[0]->isParUpdatable())
          {
            computeOneJacob(RIGHT);
          }
        }
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }

  //// Loop over the elements to calculate the divergence of the continuous flux

  // loop over element types, for the moment there should only be one
  const CFuint nbrElemTypes = elemType->size();
  cf_assert(nbrElemTypes == 1);
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // create blockaccumulator
    m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoDataCell.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // if the states in the cell are parallel updatable or the gradients need to be computed, set the cell data
      if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
      {
        // set the cell data
        setCellData();
      }

      // if the states in the cell are parallel updatable, compute the divergence of the discontinuous flx (-divFD+divhFD)
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // compute the divergence of the discontinuous flux (-divFD+divhFD)
        computeDivDiscontFlx(m_divContFlx);

        // update RHS
        updateRHS();
      }

      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
        computeGradients();
      }

      // compute the contribution to the volume numerical jacobian (Fix FJ1: gate with freeze check)
      {
        const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
        const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
        const CFuint interval = iter - iterFreeze;
        if (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0)
        {
          if ((*m_cellStates)[0]->isParUpdatable())
          {
            // add the contributions to the Jacobian (uses subcell FV blended computeDivDiscontFlx)
            computeJacobConvCorrection();
          }
        }
      }

      // release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::updateWaveSpeed()
{
  // Blended wave speed factor: (1-alpha)*(2p+1) + alpha*1
  // FR penalty factor is (2p+1), first-order FV factor is 1
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> alphaData = socket_alpha.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFreal alpha = alphaData[(*m_states[iSide])[0]->getLocalID()];
    const CFreal factor = (1.0 - alpha) * (2.0 * m_order + 1.0) + alpha * 1.0;

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();
      updateCoeff[solID] += m_waveSpeedUpd[iSide] * factor;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::computeSubcellFVResidual()
{
  // reset FV residuals
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_subcellFVRes[iSol] = 0.0;
  }

  if (!m_tensorBased || m_nbrSolPnts1D < 2)
  {
    return;
  }

  for (CFuint dir = 0; dir < m_dim; ++dir)
  {
    if (m_subcellIfaceCoords[dir].empty())
      continue;

    // compute mapped normals for all interfaces in this direction
    std::vector< RealVector > normals =
      m_cell->computeMappedCoordPlaneNormalAtMappedCoords(
        m_subcellIfaceDimList[dir], m_subcellIfaceCoords[dir]);

    const CFuint nIface = m_subcellIfaceCoords[dir].size();
    for (CFuint k = 0; k < nIface; ++k)
    {
      const CFuint idxL = m_subcellIfaceLeftIdx[dir][k];
      const CFuint idxR = m_subcellIfaceRightIdx[dir][k];

      // SFV-4: Normalize mapped normal before Riemann flux (see ConvRHSSubcellFVVS)
      const CFreal normMag = normals[k].norm2();
      m_subcellUnitNormal = normals[k] / normMag;

      RealVector& fstar = m_riemannFluxComputer->computeFlux(
        *(*m_cellStates)[idxL], *(*m_cellStates)[idxR], m_subcellUnitNormal);

      const CFreal invWL = m_subcellIfaceInvWidthL[dir][k];
      const CFreal invWR = m_subcellIfaceInvWidthR[dir][k];

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_subcellFVRes[idxL][iEq] -= fstar[iEq] * normMag * invWL;
        m_subcellFVRes[idxR][iEq] += fstar[iEq] * normMag * invWR;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  if (!m_tensorBased)
  {
    ConvRHSJacobFluxReconstruction::computeDivDiscontFlx(residuals);
    return;
  }

  DataHandle< CFreal > output = socket_alpha.getDataHandle();
  CFreal alpha = output[((*m_cellStates)[0])->getLocalID()];

  if (alpha <= m_alphaEps)
  {
    ConvRHSJacobFluxReconstruction::computeDivDiscontFlx(residuals);
    return;
  }

  // reset the extrapolated fluxes
  for (CFuint iFlxPnt = 0; iFlxPnt < m_flxPntsLocalCoords->size(); ++iFlxPnt)
  {
    m_extrapolatedFluxes[iFlxPnt] = 0.0;
  }

  // Loop over solution points to calculate the discontinuous flux.
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    m_updateVarSet->computePhysicalData(*(*m_cellStates)[iSolPnt], m_pData);

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim + m_ndimplus; ++iDim)
    {
      m_contFlx[iSolPnt][iDim] = m_updateVarSet->getFlux()(m_pData, m_cellFluxProjVects[iDim][iSolPnt]);
    }

    // extrapolate the fluxes to the flux points
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];
      const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
      m_extrapolatedFluxes[flxIdx] += (*m_solPolyValsAtFlxPnts)[flxIdx][iSolPnt] * (m_contFlx[iSolPnt][dim]);
    }
  }

  // --- FR volume term only ---
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    residuals[iSolPnt] = 0.0;
    for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
    {
      const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[iSolPnt][iDir][jSolIdx];
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          residuals[iSolPnt][iEq] -= polyCoef * (m_contFlx[jSolIdx][iDir][iEq]);
        }
      }
    }
  }

  // --- Subcell FV volume term ---
  computeSubcellFVResidual();

  // --- Blend volume term ---
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      residuals[iSolPnt][iEq] =
        (1.0 - alpha) * residuals[iSolPnt][iEq] +
        alpha * m_subcellFVRes[iSolPnt][iEq];
    }
  }

  // --- SFV-5: FR divhFD term scaled by (1-alpha) ---
  // In the Hennemann scheme, R = (1-alpha)*R_FR + alpha*R_FV.
  // The divhFD term is part of R_FR, so it must be scaled by (1-alpha).
  // The face correction in execute() is NOT scaled (see comment there).
  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];
      const CFreal divh = m_corrFctDiv[iSolPnt][flxIdx];
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        residuals[iSolPnt][iVar] -= -m_extrapolatedFluxes[flxIdx][iVar] * divh * oneMinusAlpha;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::computeJacobConvCorrection()
{
  const CFreal resFactor = getMethodData().getResFactor();
  BlockAccumulator& acc = *m_acc;

  // Get alpha for current cell
  DataHandle< CFreal > alphaData = socket_alpha.getDataHandle();
  const CFreal alpha = m_tensorBased ?
    alphaData[(*m_cellStates)[0]->getLocalID()] : 0.0;

  // Set block row and column indices
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol, (*m_cellStates)[iSol]->getLocalID());
  }

  // ==========================================================================
  // Step 1: FR Jacobian via FD, scaled by (1-alpha)
  // The parent's computeDivDiscontFlx gives R_FR = -divFD + divhFD (smooth).
  // The blended residual is (1-α)*R_FR + α*R_FV, so the FR contribution
  // to the Jacobian is (1-α)*∂R_FR/∂U.
  // ==========================================================================
  const CFreal oneMinusAlpha = 1.0 - alpha;

  // Recompute FR-only baseline
  ConvRHSFluxReconstruction::computeDivDiscontFlx(m_divContFlx);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[0][m_nbrEqs * iSol + iVar] = m_divContFlx[iSol][iVar];
    }
  }

  storeBackupsCell();

  for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
  {
    State& pertState = *(*m_cellStates)[m_pertSol];

    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      m_numJacob->perturb(m_pertVar, pertState[m_pertVar]);

      ConvRHSFluxReconstruction::computeDivDiscontFlx(m_pertDivContFlx[LEFT]);

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_pertResUpdates[0][m_nbrEqs * iSol + iVar] = m_pertDivContFlx[LEFT][iSol][iVar];
        }
      }

      m_numJacob->computeDerivative(m_pertResUpdates[0], m_resUpdates[0], m_derivResUpdates);

      // AJ-1: Scale by (1-α) instead of 1.0
      m_derivResUpdates *= (resFactor * oneMinusAlpha);

      for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
      {
        const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
        acc.addValues(jSolIdx, m_pertSol, m_pertVar, &m_derivResUpdates[m_nbrEqs * jSolIdx]);
      }

      m_numJacob->restore(pertState[m_pertVar]);

      for (CFuint iFlxPnt = 0; iFlxPnt < m_flxPntsLocalCoords->size(); ++iFlxPnt)
      {
        m_extrapolatedFluxes[iFlxPnt] = m_extrapolatedFluxesBackup[iFlxPnt];
      }
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        for (CFuint iDim = 0; iDim < m_dim + m_ndimplus; ++iDim)
        {
          m_contFlx[iSol][iDim] = m_contFlxBackup[iSol][iDim];
        }
      }
    }
  }

  // ==========================================================================
  // Step 2: Analytical subcell FV Jacobian, scaled by α
  // Adds α * ∂R_FV/∂U to the accumulator. The AUSM+ flux Jacobian is
  // computed analytically (no FD perturbation of the Riemann solver),
  // eliminating the Mach-split noise that caused the sawtooth pattern.
  // ==========================================================================
  if (m_tensorBased && alpha > m_alphaEps)
  {
    computeAnalyticalSubcellFVJacobianVS(acc, resFactor, alpha);
  }

  if (getMethodData().doComputeJacobian())
  {
    m_lss->getMatrix()->addValues(acc);
  }

  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::computeAnalyticalSubcellFVJacobianVS(
  BlockAccumulator& acc, const CFreal resFactor, const CFreal alpha)
{
  // Cache physical data and its Jacobian ∂pdata/∂U for all solution points.
  // This avoids redundant FD perturbations when the same sol point appears
  // at multiple subcell interfaces.
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_updateVarSet->computePhysicalData(*(*m_cellStates)[iSol], m_pdataPerSol[iSol]);
    computePhysicalDataJacobianVS(*(*m_cellStates)[iSol], m_pdataPerSol[iSol], m_dPdUPerSol[iSol]);
  }

  const CFreal scaleFactor = alpha * resFactor;

  // Loop over subcell interfaces in each direction
  for (CFuint dir = 0; dir < m_dim; ++dir)
  {
    if (m_subcellIfaceCoords[dir].empty()) continue;

    // Compute mapped normals for all interfaces in this direction
    std::vector< RealVector > normals =
      m_cell->computeMappedCoordPlaneNormalAtMappedCoords(
        m_subcellIfaceDimList[dir], m_subcellIfaceCoords[dir]);

    const CFuint nIface = m_subcellIfaceCoords[dir].size();
    for (CFuint k = 0; k < nIface; ++k)
    {
      const CFuint idxL = m_subcellIfaceLeftIdx[dir][k];
      const CFuint idxR = m_subcellIfaceRightIdx[dir][k];

      const CFreal normMag = normals[k].norm2();
      m_subcellUnitNormal = normals[k] / normMag;

      const CFreal invWL = m_subcellIfaceInvWidthL[dir][k];
      const CFreal invWR = m_subcellIfaceInvWidthR[dir][k];

      // Compute analytical AUSM+ flux Jacobian ∂F*/∂U_L and ∂F*/∂U_R
      computeAUSMPlusFluxJacobianVS(
        m_pdataPerSol[idxL], m_pdataPerSol[idxR], m_subcellUnitNormal,
        m_dPdUPerSol[idxL], m_dPdUPerSol[idxR],
        m_dFluxdUL, m_dFluxdUR);

      // Accumulate into Jacobian matrix:
      //   R_FV[L] -= F* · normMag · invWL
      //   R_FV[R] += F* · normMag · invWR
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        for (CFuint jVar = 0; jVar < m_nbrEqs; ++jVar)
        {
          const CFreal jL = m_dFluxdUL(iEq, jVar) * normMag;
          const CFreal jR = m_dFluxdUR(iEq, jVar) * normMag;

          // ∂R_FV[L]/∂U_L, ∂R_FV[L]/∂U_R  (minus sign from R_FV[L] -= ...)
          acc.addValue(idxL, idxL, iEq, jVar, -scaleFactor * jL * invWL);
          acc.addValue(idxL, idxR, iEq, jVar, -scaleFactor * jR * invWL);

          // ∂R_FV[R]/∂U_L, ∂R_FV[R]/∂U_R  (plus sign from R_FV[R] += ...)
          acc.addValue(idxR, idxL, iEq, jVar,  scaleFactor * jL * invWR);
          acc.addValue(idxR, idxR, iEq, jVar,  scaleFactor * jR * invWR);
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::computePhysicalDataJacobianVS(
  State& state, const RealVector& pdataBase, RealMatrix& dPdU)
{
  // FD perturbation of computePhysicalData to get ∂pdata/∂U.
  // This is noise-free because computePhysicalData is a smooth function
  // (thermodynamics, algebra — no sign/max/min branches).
  for (CFuint jVar = 0; jVar < m_nbrEqs; ++jVar)
  {
    const CFreal origVal = state[jVar];
    const CFreal eps = m_pdataFDEps * std::max(std::abs(origVal), 1.0);

    state[jVar] = origVal + eps;
    m_updateVarSet->computePhysicalData(state, m_pdataPerturbed);
    state[jVar] = origVal;

    const CFreal invEps = 1.0 / eps;
    for (CFuint p = 0; p < m_pdataSize; ++p)
    {
      dPdU(p, jVar) = (m_pdataPerturbed[p] - pdataBase[p]) * invEps;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobSubcellFVVS::computeAUSMPlusFluxJacobianVS(
  const RealVector& pdataL, const RealVector& pdataR,
  const RealVector& n,
  const RealMatrix& dPdUL, const RealMatrix& dPdUR,
  RealMatrix& dFluxdUL, RealMatrix& dFluxdUR)
{
  // EulerTerm physical data indices
  const CFuint iRHO = 0, iP = 1, iH = 2, iA = 4, iVX = 7;

  // Extract physical quantities
  const CFreal rhoL = pdataL[iRHO], rhoR = pdataR[iRHO];
  const CFreal pL   = pdataL[iP],   pR   = pdataR[iP];
  const CFreal HL   = pdataL[iH],   HR   = pdataR[iH];
  const CFreal aL   = pdataL[iA],   aR   = pdataR[iA];

  // Normal velocities
  CFreal unL = 0.0, unR = 0.0;
  for (CFuint d = 0; d < m_dim; ++d)
  {
    unL += pdataL[iVX + d] * n[d];
    unR += pdataR[iVX + d] * n[d];
  }

  // Interface sound speed (choiceA12 = 3: geometric mean)
  const CFreal a12 = std::sqrt(aL * aR);
  cf_assert(a12 > 0.0);

  // Mach numbers
  const CFreal ML = unL / a12;
  const CFreal MR = unR / a12;

  // AUSM+ parameters (standard defaults)
  const CFreal beta = 1.0 / 8.0;
  const CFreal alpha_coeff = 3.0 / 16.0;

  // ---- Mach number splitting M4± and derivatives ----
  CFreal M4p, M4m, dM4pdM, dM4mdM;

  if (std::abs(ML) >= 1.0) {
    M4p = 0.5 * (ML + std::abs(ML));
    dM4pdM = (ML > 0.0) ? 1.0 : 0.0;
  } else {
    const CFreal MLp1 = ML + 1.0;
    const CFreal ML2m1 = ML * ML - 1.0;
    M4p = 0.25 * MLp1 * MLp1 + beta * ML2m1 * ML2m1;
    dM4pdM = 0.5 * MLp1 + 4.0 * beta * ML * ML2m1;
  }

  if (std::abs(MR) >= 1.0) {
    M4m = 0.5 * (MR - std::abs(MR));
    dM4mdM = (MR < 0.0) ? 1.0 : 0.0;
  } else {
    const CFreal MRm1 = MR - 1.0;
    const CFreal MR2m1 = MR * MR - 1.0;
    M4m = -0.25 * MRm1 * MRm1 - beta * MR2m1 * MR2m1;
    dM4mdM = -0.5 * MRm1 - 4.0 * beta * MR * MR2m1;
  }

  const CFreal m12 = M4p + M4m;
  const bool upwindL = (m12 > 0.0);
  const CFreal rhoUp = upwindL ? rhoL : rhoR;
  const CFreal mflux = a12 * m12 * rhoUp;

  // ---- Pressure splitting P5± and derivatives ----
  CFreal P5p, P5m, dP5pdM, dP5mdM;

  if (std::abs(ML) >= 1.0) {
    P5p = (ML > 0.0) ? 1.0 : 0.0;
    dP5pdM = 0.0;
  } else {
    const CFreal MLp1 = ML + 1.0;
    const CFreal ML2m1 = ML * ML - 1.0;
    P5p = 0.25 * MLp1 * MLp1 * (2.0 - ML)
        + alpha_coeff * ML * ML2m1 * ML2m1;
    dP5pdM = 0.5 * MLp1 * (2.0 - ML) - 0.25 * MLp1 * MLp1
           + alpha_coeff * ML2m1 * (5.0 * ML * ML - 1.0);
  }

  if (std::abs(MR) >= 1.0) {
    P5m = (MR < 0.0) ? 1.0 : 0.0;
    dP5mdM = 0.0;
  } else {
    const CFreal MRm1 = MR - 1.0;
    const CFreal MR2m1 = MR * MR - 1.0;
    P5m = 0.25 * MRm1 * MRm1 * (2.0 + MR)
        - alpha_coeff * MR * MR2m1 * MR2m1;
    dP5mdM = 0.5 * MRm1 * (2.0 + MR) + 0.25 * MRm1 * MRm1
           - alpha_coeff * MR2m1 * (5.0 * MR * MR - 1.0);
  }

  // ---- Intermediate derivatives ----
  // ∂a12/∂aL = 0.5*a12/aL, ∂a12/∂aR = 0.5*a12/aR
  const CFreal da12_daL = 0.5 * a12 / aL;
  const CFreal da12_daR = 0.5 * a12 / aR;

  // ∂ML/∂unL = 1/a12
  // ∂ML/∂aL = -ML * da12_daL / a12 = -0.5*ML/aL  (through a12)
  // ∂ML/∂aR = -0.5*ML/aR
  const CFreal dML_dunL = 1.0 / a12;
  const CFreal dML_daL  = -0.5 * ML / aL;
  const CFreal dML_daR  = -0.5 * ML / aR;
  const CFreal dMR_dunR = 1.0 / a12;
  const CFreal dMR_daL  = -0.5 * MR / aL;
  const CFreal dMR_daR  = -0.5 * MR / aR;

  // ∂m12/∂(L quantities): through M4+(ML)
  // ∂m12/∂(R quantities): through M4-(MR)
  // Also both M4+ and M4- depend on a12 (through ML, MR)
  const CFreal dm12_dunL = dM4pdM * dML_dunL;
  const CFreal dm12_dunR = dM4mdM * dMR_dunR;
  const CFreal dm12_daL  = dM4pdM * dML_daL + dM4mdM * dMR_daL;
  const CFreal dm12_daR  = dM4pdM * dML_daR + dM4mdM * dMR_daR;

  // ∂mflux/∂(...) = ∂(a12*m12*rhoUp)/∂(...)
  const CFreal dmflux_dunL = dM4pdM * rhoUp;  // a12 cancels: a12 * dM4p/dM * 1/a12
  const CFreal dmflux_dunR = dM4mdM * rhoUp;
  const CFreal dmflux_daL  = (da12_daL * m12 + a12 * dm12_daL) * rhoUp;
  const CFreal dmflux_daR  = (da12_daR * m12 + a12 * dm12_daR) * rhoUp;
  const CFreal dmflux_drhoL = upwindL ? (a12 * m12) : 0.0;
  const CFreal dmflux_drhoR = upwindL ? 0.0 : (a12 * m12);

  // ∂p12/∂(...)  where p12 = P5p*pL + P5m*pR
  const CFreal dp12_dpL   = P5p;
  const CFreal dp12_dpR   = P5m;
  const CFreal dp12_dunL  = dP5pdM * dML_dunL * pL;
  const CFreal dp12_dunR  = dP5mdM * dMR_dunR * pR;
  const CFreal dp12_daL   = dP5pdM * dML_daL * pL + dP5mdM * dMR_daL * pR;
  const CFreal dp12_daR   = dP5pdM * dML_daR * pL + dP5mdM * dMR_daR * pR;

  // Upwind quantities
  const CFreal HUp = upwindL ? HL : HR;

  // ====================================================================
  // Build ∂F/∂U_L and ∂F/∂U_R via chain rule:
  //   ∂F/∂U[jVar] = Σ_q ∂F/∂pdata[q] × ∂pdata[q]/∂U[jVar]
  // Only the entries for q ∈ {RHO, VX..VZ, P, A, H, firstScalar+species}
  // are non-zero, so the sum is sparse.
  // ====================================================================
  dFluxdUL = 0.0;
  dFluxdUR = 0.0;

  for (CFuint jVar = 0; jVar < m_nbrEqs; ++jVar)
  {
    // --- LEFT side chain rule components ---
    CFreal dunL_dUj = 0.0;
    for (CFuint d = 0; d < m_dim; ++d)
      dunL_dUj += n[d] * dPdUL(iVX + d, jVar);

    const CFreal daL_dUj   = dPdUL(iA, jVar);
    const CFreal drhoL_dUj = dPdUL(iRHO, jVar);
    const CFreal dpL_dUj   = dPdUL(iP, jVar);
    const CFreal dHL_dUj   = dPdUL(iH, jVar);

    // ∂mflux/∂U_L[j]
    const CFreal dmflux_dUjL = dmflux_dunL * dunL_dUj
                             + dmflux_daL  * daL_dUj
                             + dmflux_drhoL * drhoL_dUj;
    // ∂p12/∂U_L[j]
    const CFreal dp12_dUjL = dp12_dunL * dunL_dUj
                           + dp12_daL  * daL_dUj
                           + dp12_dpL  * dpL_dUj;

    // --- RIGHT side chain rule components ---
    CFreal dunR_dUj = 0.0;
    for (CFuint d = 0; d < m_dim; ++d)
      dunR_dUj += n[d] * dPdUR(iVX + d, jVar);

    const CFreal daR_dUj   = dPdUR(iA, jVar);
    const CFreal drhoR_dUj = dPdUR(iRHO, jVar);
    const CFreal dpR_dUj   = dPdUR(iP, jVar);
    const CFreal dHR_dUj   = dPdUR(iH, jVar);

    const CFreal dmflux_dUjR = dmflux_dunR * dunR_dUj
                             + dmflux_daR  * daR_dUj
                             + dmflux_drhoR * drhoR_dUj;
    const CFreal dp12_dUjR = dp12_dunR * dunR_dUj
                           + dp12_daR  * daR_dUj
                           + dp12_dpR  * dpR_dUj;

    // --- Species/continuity equations ---
    for (CFuint iSpec = 0; iSpec < m_nbSpecies; ++iSpec)
    {
      const CFuint iEq = iSpec;
      const CFuint iScalar = m_firstScalarPdata + iSpec;

      if (iScalar < m_pdataSize)
      {
        // Multi-species (CNEQ/TCNEQ): F[iSpec] = mflux * Y_iSpec_upwind
        const CFreal YiUp = upwindL ? pdataL[iScalar] : pdataR[iScalar];

        dFluxdUL(iEq, jVar) = dmflux_dUjL * YiUp;
        if (upwindL) dFluxdUL(iEq, jVar) += mflux * dPdUL(iScalar, jVar);

        dFluxdUR(iEq, jVar) = dmflux_dUjR * YiUp;
        if (!upwindL) dFluxdUR(iEq, jVar) += mflux * dPdUR(iScalar, jVar);
      }
      else
      {
        // Single-species (Euler): F[0] = mflux (continuity, Y=1)
        dFluxdUL(iEq, jVar) = dmflux_dUjL;
        dFluxdUR(iEq, jVar) = dmflux_dUjR;
      }
    }

    // --- Momentum equations: F[mom_d] = mflux * vel_d_upwind + p12 * n[d] ---
    for (CFuint d = 0; d < m_dim; ++d)
    {
      const CFuint iEq = m_eulerStartIdx + d;
      const CFreal velUp = upwindL ? pdataL[iVX + d] : pdataR[iVX + d];

      dFluxdUL(iEq, jVar) = dmflux_dUjL * velUp + dp12_dUjL * n[d];
      if (upwindL) dFluxdUL(iEq, jVar) += mflux * dPdUL(iVX + d, jVar);

      dFluxdUR(iEq, jVar) = dmflux_dUjR * velUp + dp12_dUjR * n[d];
      if (!upwindL) dFluxdUR(iEq, jVar) += mflux * dPdUR(iVX + d, jVar);
    }

    // --- Energy equation: F[energy] = mflux * H_upwind ---
    {
      const CFuint iEq = m_eulerStartIdx + m_dim;

      dFluxdUL(iEq, jVar) = dmflux_dUjL * HUp;
      if (upwindL) dFluxdUL(iEq, jVar) += mflux * dHL_dUj;

      dFluxdUR(iEq, jVar) = dmflux_dUjR * HUp;
      if (!upwindL) dFluxdUR(iEq, jVar) += mflux * dHR_dUj;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
