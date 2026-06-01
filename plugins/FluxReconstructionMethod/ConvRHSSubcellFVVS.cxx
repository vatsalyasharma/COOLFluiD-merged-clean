// Copyright (C) 2026
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMethod/ConvRHSSubcellFVVS.hh"
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

MethodCommandProvider< ConvRHSSubcellFVVS,
                       FluxReconstructionSolverData,
                       FluxReconstructionModule >
  ConvRHSSubcellFVVSProvider("ConvRHSSubcellFVVS");

//////////////////////////////////////////////////////////////////////////////

void ConvRHSSubcellFVVS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("AlphaEps",
    "Skip subcell FV when alpha <= AlphaEps (default 1e-12).");
}

//////////////////////////////////////////////////////////////////////////////

ConvRHSSubcellFVVS::ConvRHSSubcellFVVS(const std::string& name) :
  ConvRHSFluxReconstruction(name),
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

void ConvRHSSubcellFVVS::configure ( Config::ConfigArgs& args )
{
  ConvRHSFluxReconstruction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
ConvRHSSubcellFVVS::needsSockets()
{
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > result =
    ConvRHSFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSSubcellFVVS::setup()
{
  CFAUTOTRACE;

  ConvRHSFluxReconstruction::setup();

  // resize subcell FV residual storage
  m_subcellFVRes.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_subcellFVRes[iSol].resize(m_nbrEqs);
  }

  // unit normal for subcell FV Riemann flux
  m_subcellUnitNormal.resize(m_dim);

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
    CFLog(WARN, "ConvRHSSubcellFVVS: non tensor-product element ("
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

void ConvRHSSubcellFVVS::unsetup()
{
  ConvRHSFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSSubcellFVVS::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "ConvRHSSubcellFVVS::execute()\n");

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

      // release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSSubcellFVVS::updateWaveSpeed()
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

void ConvRHSSubcellFVVS::computeSubcellFVResidual()
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

      // SFV-4: The mapped normals are cofactor vectors with magnitude = face area element.
      // The Riemann solver (AUSM+, LaxFriedrichs, etc.) expects a UNIT normal because
      // its Mach splitting M=u·n/a is nonlinear: M+(|n|*M_phys) != |n|*M+(M_phys).
      // With raw mapped normals (|n|~h/2), hypersonic flows appear subsonic to AUSM.
      // Fix: normalize before Riemann, scale result by |n| after.
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

void ConvRHSSubcellFVVS::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  if (!m_tensorBased)
  {
    ConvRHSFluxReconstruction::computeDivDiscontFlx(residuals);
    return;
  }

  DataHandle< CFreal > output = socket_alpha.getDataHandle();
  CFreal alpha = output[((*m_cellStates)[0])->getLocalID()];

  if (alpha <= m_alphaEps)
  {
    ConvRHSFluxReconstruction::computeDivDiscontFlx(residuals);
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

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

