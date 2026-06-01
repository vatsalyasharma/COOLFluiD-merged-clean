// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// SFV-6B: Boundary diffusive corrections — NO alpha scaling.
// Per Rueda-Ramírez et al. (2021), diffusion is discretized entirely with
// high-order FR (no FV blending). Volume diffusion (DiffRHSBlendingVS) already
// has no scaling. Boundary corrections must match to avoid an inconsistent
// operator (full volume diffusion + reduced boundary diffusion).

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSFluxReconstructionBlendingVS.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSFluxReconstructionBlendingVS,
                       FluxReconstructionSolverData,
                       FluxReconstructionModule >
DiffBndCorrectionsRHSFluxReconstructionBlendingVSProvider("DiffBndCorrectionsRHSBlendingVS");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSFluxReconstructionBlendingVS::
DiffBndCorrectionsRHSFluxReconstructionBlendingVS(const std::string& name) :
  DiffBndCorrectionsRHSFluxReconstruction(name),
  socket_alpha("alpha")
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
DiffBndCorrectionsRHSFluxReconstructionBlendingVS::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result =
    DiffBndCorrectionsRHSFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstructionBlendingVS::updateRHS()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // SFV-6B: No alpha scaling — full FR diffusion at boundaries,
  // matching the unblended volume diffusion (SFV-6).
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      rhs[resID+iVar] += resFactor*m_corrections[iState][iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
