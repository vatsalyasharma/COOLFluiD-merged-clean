// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"


#include "StdUpdateSol.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Common/BadValueException.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/ConvectiveVarSet.hh"

#include "Framework/ConvergenceMethod.hh"
#include "Framework/ConvergenceMethodData.hh"
#include "Framework/ConvergenceStatus.hh"
#include "NewtonMethod/NewtonIteratorData.hh"
#include "NewtonMethod/NewtonIterator.hh"

#include "Framework/PathAppender.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/ElementTypeData.hh"
#include "Framework/StdTrsGeoBuilder.hh"

#include <algorithm>
#include <cmath>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUpdateSol, NewtonIteratorData, NewtonMethodModule> 
stdUpdateSolProvider("StdUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void StdUpdateSol::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption<vector<CFreal>,Config::DynamicOption<> >("Relaxation","Relaxation factor");
  options.addConfigOption< bool >("Validate","Check that each update creates variables with physical meaning");
  options.addConfigOption< bool >("PositivityLineSearch","Enable coupled positivity-preserving line search on Newton update");
  options.addConfigOption<vector<CFreal> >("PositivityMinValues","Minimum values per equation for positivity line search");
  options.addConfigOption< bool,Config::DynamicOption<> >("FilterNewtonUpdate","Damp Newton update modes using frozen entropy-filter strength");
}

//////////////////////////////////////////////////////////////////////////////

StdUpdateSol::StdUpdateSol(const std::string& name) : 
  NewtonIteratorCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_invalidStates("invalidStates"),
  socket_filterZeta("outputPP")
{
  addConfigOptionsTo(this);
   
  m_alpha = vector<CFreal>();
  setParameter("Relaxation",&m_alpha);

  m_validate = false;
  setParameter("Validate",&m_validate);

  m_positivityLineSearch = false;
  setParameter("PositivityLineSearch",&m_positivityLineSearch);

  m_positivityMinValues = vector<CFreal>();
  setParameter("PositivityMinValues",&m_positivityMinValues);

  m_filterNewtonUpdate = false;
  setParameter("FilterNewtonUpdate",&m_filterNewtonUpdate);
}

//////////////////////////////////////////////////////////////////////////////

void StdUpdateSol::setup()
{
  const CFuint nbStates = socket_states.getDataHandle().size();

  if (m_validate) {
    socket_invalidStates.getDataHandle().resize(nbStates);
  }
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  if (m_alpha.size() == 0) {
    m_alpha.resize(nbEqs);
    for (CFuint i = 0; i < nbEqs; ++i) {
      m_alpha[i] = 1.;
    }
  }

  if (m_alpha.size() == 1) {
    const CFreal value = m_alpha[0];
    m_alpha.resize(nbEqs);
    for (CFuint i = 0; i < nbEqs; ++i) {
      m_alpha[i] = value;
    }
  }

  if (m_alpha.size() != nbEqs) {
    throw BadValueException (FromHere(),"StdUpdateSol::setup() : m_alpha.size() != nbEqs");
  }

  if (m_positivityLineSearch) {
    if (m_positivityMinValues.size() != nbEqs) {
      throw BadValueException (FromHere(),
        "StdUpdateSol::setup() : PositivityMinValues size must equal nbEqs when PositivityLineSearch is enabled");
    }

    CFLog(INFO, "StdUpdateSol: positivity line search enabled with minValues = [");
    for (CFuint i = 0; i < nbEqs; ++i) {
      CFLog(INFO, " " << m_positivityMinValues[i]);
    }
    CFLog(INFO, " ]\n");
  }

  if (m_filterNewtonUpdate) {
    m_geoBuilder.setup();
    CFLog(INFO, "StdUpdateSol: filtered Newton update enabled\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdUpdateSol::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  // rhs is the temporary placeholder for the dU
  DataHandle<CFreal>& dU = rhs;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint states_size = states.size();

  SafePtr<FilterState> filterState = getMethodData().getFilterState();
  SafePtr<FilterRHS> filterRHS     = getMethodData().getFilterRHS();

  if (m_filterNewtonUpdate)
  {
    DataHandle<CFreal> filterZeta = socket_filterZeta.getDataHandle();

    SafePtr<TopologicalRegionSet> cells =
      MeshDataStack::getActive()->getTrs("InnerCells");
    SafePtr<vector<ElementTypeData> > elemType =
      MeshDataStack::getActive()->getElementTypeData();

    StdTrsGeoBuilder::GeoData& geoData =
      m_geoBuilder.getDataGE();
    geoData.trs = cells;

    CFuint nCellsFiltered = 0;
    const CFuint nbrElemTypes = elemType->size();

    for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
    {
      const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
      const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

      for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
      {
        geoData.idx = elemIdx;
        GeometricEntity* cell = m_geoBuilder.buildGE();
        vector<State*>* cellStates = cell->getStates();
        const CFuint nbCellStates = cellStates->size();

        CFreal zetaCell = 0.0;
        for (CFuint iSol = 0; iSol < nbCellStates; ++iSol)
        {
          const CFuint stateID = (*((*cellStates)[iSol])).getLocalID();
          if (stateID < filterZeta.size())
          {
            const CFreal zeta = filterZeta[stateID];
            if (zeta > zetaCell) zetaCell = zeta;
          }
        }

        if (zetaCell > 1.0e-10)
        {
          ++nCellsFiltered;

          // Vatsalya: damp non-mean Newton-update modes using frozen entropy-filter strength.
          const CFreal sigma = std::exp(-zetaCell);
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
          {
            CFreal meanDU = 0.0;
            for (CFuint iSol = 0; iSol < nbCellStates; ++iSol)
            {
              const CFuint stateID = (*((*cellStates)[iSol])).getLocalID();
              meanDU += dU(stateID, iEq, nbEqs);
            }
            meanDU /= static_cast<CFreal>(nbCellStates);

            for (CFuint iSol = 0; iSol < nbCellStates; ++iSol)
            {
              const CFuint stateID = (*((*cellStates)[iSol])).getLocalID();
              const CFreal orig = dU(stateID, iEq, nbEqs);
              dU(stateID, iEq, nbEqs) = meanDU + sigma * (orig - meanDU);
            }
          }
        }

        m_geoBuilder.releaseGE();
      }
    }

    if (nCellsFiltered > 0) {
      CFLog(INFO, "StdUpdateSol: filtered Newton update damped dU in "
            << nCellsFiltered << " cells\n");
    }
  }

  CFuint nClipped = 0;
  
  for (CFuint iState = 0; iState < states_size; ++iState)
  {
    State& cur_state = *states[iState];
    // do the update only if the state is parallel updatable
    if (cur_state.isParUpdatable())
    {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
	filterRHS->filter(iEq, dU(iState, iEq, nbEqs));
      }

      CFreal omega = 1.0;
      if (m_positivityLineSearch)
      {
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        {
          const CFreal minVal = m_positivityMinValues[iEq];
          if (minVal <= -1.0e29) continue;

          const CFreal step = m_alpha[iEq] * dU(iState, iEq, nbEqs);
          const CFreal newVal = cur_state[iEq] + step;
          if (newVal < minVal && step < 0.0)
          {
            const CFreal omegaEq = 0.99 * (minVal - cur_state[iEq]) / step;
            omega = std::min(omega, std::max(omegaEq, 0.0));
          }
        }
        if (omega < 1.0) ++nClipped;
      }

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
	cur_state[iEq] += omega * m_alpha[iEq] * dU(iState, iEq, nbEqs);
      }
      
      // apply a polymorphic filter to the state
      filterState->filter(cur_state);
    }
    else {
      // reset to 0 the RHS for ghost states in order to avoid 
      // inconsistencies in the parallel L2 norm computation
      cf_assert(!cur_state.isParUpdatable());
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
	dU(iState, iEq, nbEqs) = 0.;
      }
    }
  }

  if (m_positivityLineSearch && nClipped > 0) {
    CFLog(INFO, "StdUpdateSol: positivity line search clipped " << nClipped << " states\n");
  }
  
  ///This vector will temporally contain the ID of any vector considered not valid.
  std::vector<CFuint> badStatesIDs;

  ///Name of the file where to store invalid states in.
  boost::filesystem::path filepath ( "unphysical_states.plt");

  filepath = PathAppender::getInstance().appendParallel( filepath );

  // loop for unphysicalness check:
  if ( m_validate )
  {
    Common::SafePtr<SpaceMethod> theSpaceMethod = getMethodData().getCollaborator<SpaceMethod>();
    Common::SafePtr<SpaceMethodData> theSpaceMethodData = theSpaceMethod->getSpaceMethodData();
    Common::SafePtr<Framework::ConvectiveVarSet> theVarSet = theSpaceMethodData->getUpdateVar();

    ///Reset the vector:
    badStatesIDs.resize(0);

    for (CFuint iState = 0; iState < states_size; ++iState)
    {
      State& cur_state = *states[iState];
      if ( !( theVarSet->isValid(cur_state) ) )
      {
        badStatesIDs.push_back( iState );
      }
    } // for all states
  } // if validate

  if (badStatesIDs.size() > 0)
  {
    cf_assert(m_validate);
    // sort the invalid states in order to be able to search later on 
    sort(badStatesIDs.begin(), badStatesIDs.end());
    
    DataHandle<CFreal> invalidStates = socket_invalidStates.getDataHandle();
    cf_assert(invalidStates.size() == states_size);
    invalidStates = 0.;
    for (CFuint i = 0; i < badStatesIDs.size(); ++i) {
      invalidStates[badStatesIDs[i]] = 1.;
    }
    
    CFLog(VERBOSE, "StdUpdateSol::execute() => [" << badStatesIDs.size() << "] invalid states detected\n");
    // correct all unphysical states 
    correctUnphysicalStates(badStatesIDs);
    
    ///Gets the global iteration
    CFuint cur_global_iter = SubSystemStatusStack::getActive()->getNbIter();
    ///Gets the Newton method procedure iteration
    CFuint cur_iter = getMethodData().getConvergenceStatus().iter;

    ///opens a file to put the states which have unphysical variables
    ofstream fout ( filepath.string().c_str(), ios::app );

    if ( !fout.is_open() )
    {
      cout << endl << "File [" << filepath.string().c_str() << "]  is not open.\n";
    }
    ///
    const std::vector<std::string>& theVarNames=getMethodData().getCollaborator<SpaceMethod>()->getSpaceMethodData()->getUpdateVar()->getVarNames();
    CFuint dim = PhysicalModelStack::getActive()->getDim();
    ///

    ///Zone header:
    fout << "VARIABLES = \"X\"\n ";
    fout << "\"Y\"" << endl;
    if (dim == 3 ){
        fout << "\"Z\"" << endl;
    }
    for (CFuint i = 0; i < theVarNames.size(); i++){
      fout << "\"" << theVarNames[i] << "\"" << endl;
    }
    ///


    CFuint badStates_size = badStatesIDs.size();

    ///Line  
    fout << "ZONE T = \"Global iteration "<< cur_global_iter << ", Newton iteration " << cur_iter << "\"" << endl;
    ///Line
    fout << "I="<< badStates_size << ", J=1, K=1, ZONETYPE=Ordered"<< endl;
    ///Line
    fout << "DATAPACKING=POINT" << endl;
    ///Line
    fout << "DT=(";
    for ( CFuint i = 0 ; i < dim + theVarNames.size(); i++){
      fout << "DOUBLE ";
    }

    fout << ")" << endl;
    ///

    for (CFuint index = 0; index < badStates_size; ++index)
    {
     State& invalid_state  = *states[ badStatesIDs[index] ];
     fout << invalid_state.getCoordinates() << " " << invalid_state << std::endl;
    } // for all invalid states

    fout.close();

  } // writing badStates

  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdUpdateSol::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);

  if (m_filterNewtonUpdate) {
    result.push_back(&socket_filterZeta);
  }
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSource> > StdUpdateSol::providesSockets()
{
  vector<SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_invalidStates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
