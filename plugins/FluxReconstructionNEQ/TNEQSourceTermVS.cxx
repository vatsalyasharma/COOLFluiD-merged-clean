#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "Framework/PhysicalChemicalLibrary.hh"

#include "FluxReconstructionNEQ/FluxReconstructionNEQ.hh"
#include "FluxReconstructionNEQ/TNEQSourceTermVS.hh"
#include "NEQ/NEQReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQ;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TNEQSourceTermVS, FluxReconstructionSolverData, FluxReconstructionNEQModule>
TNEQSourceTermVSProvider("TNEQSourceTermVS");

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTermVS::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

TNEQSourceTermVS::TNEQSourceTermVS(const std::string& name) :
    CNEQSourceTermVS(name),
    m_omegaRad(),
    m_divV(),
    m_pe(),
    m_omegaTv(),
    m_refData(CFNULL)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

TNEQSourceTermVS::~TNEQSourceTermVS()
{
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTermVS::getSourceTermData()
{
  CNEQSourceTermVS::getSourceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTermVS::addSourceTerm(RealVector& resUpdates)
{
//   // get the datahandle of the rhs
//   DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
//
//   // get residual factor
//   const CFreal resFactor = getMethodData().getResFactor();
//
//   // loop over solution points in this cell to add the source term
//   CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrSol = m_cellStates->size();

  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();

  const CFuint iEqSS = eqSS.getEqSS();
  SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel();
  const CFuint nbSpecies = term->getNbScalarVars(0);
  const CFuint nbEvEqs = term->getNbScalarVars(1);
  const CFuint nbEulerEq = m_dim + 2;
  const CFuint nbEqs = eqSS.getNbEqsSS();
  const vector<CFuint>& varIDs = MultiScalarVarSet<Euler2DVarSet>::EULERSET::getEqSetData()[0].getEqSetVarIDs();

  bool doComputeST = false;
  if (varIDs[0] > 0 && (iEqSS == 0 && nbEqs >= nbSpecies)) {
    doComputeST = true;
  }

  if ((varIDs[0] == 0 && (iEqSS == 0) && (nbEqs >= nbEulerEq+nbSpecies)) ||
      (varIDs[0] == 0 && (iEqSS == 1)))	 {
    doComputeST = true;
  }

  if (doComputeST) {
//     // this source term is for axisymmetric flows
//     const vector<State*>* const states = element->getStates();
//
//     cf_assert(states->size() == 1);
    for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
    {
      m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);

//     // this cannot be used as is in weakly coupled simulation
//     if (_includeAxiNS) {
//       computeAxiNS(element, source);
//     }

      RealVector& refData = m_eulerVarSet->getModel()->getReferencePhysicalData();

      SafePtr<NEQReactionTerm> rt = PhysicalModelStack::getActive()->getImplementor()->
        getSourceTerm().d_castTo<Physics::NEQ::NEQReactionTerm>();

      CFreal pdim = (m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::P] + m_eulerVarSet->getModel()->getPressInf())*
        refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::P];
      cf_assert(pdim > 0.);
      CFreal Tdim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
      cf_assert(Tdim > 0.);
      CFreal rhodim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO];
      cf_assert(rhodim > 0.);

      const CFuint firstSpecies = term->getFirstScalarVar(0);
      for (CFuint i = 0; i < nbSpecies; ++i)
      {
        m_ys[i] = m_solPhysData[firstSpecies + i];
      }

      State *const currState = (*m_cellStates)[iSol];
      setVibTemperature(m_solPhysData, *currState, m_tvDim);
      m_tvDim *= refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];

      cf_assert(m_tvDim > 0.0);

//     CFLog(DEBUG_MAX, "ChemNEQST::computeSource() => T = " << Tdim << ", p = " << pdim
// 	  << ", rho = " << rhodim << ", Tv = " << _tvDim
// 	  << ", ys = [" << _ys << "], ys.sum() = " << _ys.sum() << "\n");

      m_omegaTv = 0.0;
      m_omegaRad = 0.0;

      RealMatrix jacob(m_nbrEqs,m_nbrEqs);

     // compute the conservation equation source term
     // AM: ugly but effective
     // the real solution would be to implement the function
     // MutationLibrary2OLD::getSource()
     if (this-> m_library->getName() != "Mutation2OLD"
	 && this-> m_library->getName() != "MutationPanesi"
	 && this-> m_library->getName() != "Mutationpp") {
       this-> m_library->getSource(Tdim, this-> m_tvDim, pdim, rhodim, this-> m_ys,
				  false, this-> m_omega, m_omegaTv, m_omegaRad, jacob);
      }
      else {
        // compute the mass production/destruction term
        m_library->getMassProductionTerm(Tdim, this-> m_tvDim, pdim, rhodim, this-> m_ys,
					   false, this-> m_omega, jacob);

        // compute energy relaxation and excitation term
        if (nbEvEqs > 0) {
	  // this can include all source terms for the electron energy equation if there is no vibration
	  m_library->getSourceTermVT(Tdim, this-> m_tvDim, pdim, rhodim, m_omegaTv, m_omegaRad);
        }
      }

      cf_assert(m_ys.sum() > 0.99 && m_ys.sum() < 1.0001);

    //Vatsalya: original had unconditional getMassProductionTerm() call here AFTER the if/else block that already computed it; doubled chemistry cost per cell per sol point with identical results

      CFLog(DEBUG_MAX, "ChemNEQST::computeSource() => omega = " << m_omega << "\n");

      const vector<CFuint>& speciesVarIDs = MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].getEqSetVarIDs();
      const vector<CFuint>& evVarIDs = MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[1].getEqSetVarIDs();

    //     const CFreal ovOmegaRef = PhysicalModelStack::getActive()->
    //       getImplementor()->getRefLength()/(refData[UPDATEVAR::PTERM::V]*
    // 					sourceRefData[NEQReactionTerm::TAU]);

      const CFreal ovOmegaRef = PhysicalModelStack::getActive()->getImplementor()->
        getRefLength()/(refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::V]);

      const CFreal ovOmegavRef = PhysicalModelStack::getActive()->getImplementor()->
        getRefLength()/((*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::H]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::V]);

      for (CFuint i = 0; i < nbSpecies; ++i)
      {
//         m_srcTerm[speciesVarIDs[i]] = m_omega[i]*ovOmegaRef;
	resUpdates[m_nbrEqs*iSol + speciesVarIDs[i]] = m_omega[i]*ovOmegaRef;
      }

      SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel();
      const CFuint nbSpecies = term->getNbScalarVars(0);
      const CFuint nbEvEqs = term->getNbScalarVars(1);
      const CFuint TID = nbSpecies + m_dim;
      const CFuint TED = nbSpecies + m_dim + nbEvEqs;

      if (nbEvEqs > 0) { cf_always_assert(TID == (evVarIDs[0]-1)); }

//       m_srcTerm[TID] = -m_omegaRad*ovOmegavRef;
      resUpdates[m_nbrEqs*iSol + TID] = -m_omegaRad*ovOmegavRef;

      // NOTE (Finding 8): Radiation coupling and pe*div(v) terms are NOT implemented.
      // Both RHS and Jacobian consistently omit them, so there is no mismatch.
      // To enable radiation coupling, uncomment below AND add matching Jacobian terms
      // in getSToStateJacobian().
//     if (m_hasRadiationCoupling) {
//       cf_assert(elemID < this->_qrad.size());
//       const CFreal qRad = 1.0*this->m_qrad[elemID]*ovOmegavRef;
//       m_srcTerm[TID] = - qRad;
//       if (m_library->presenceElectron()) {
//         m_srcTerm[TED] -= qRad;
//       }
//     }

      for (CFuint i = 0; i < nbEvEqs; ++i) {
//         m_srcTerm[evVarIDs[i]] = m_omegaTv[i]*ovOmegavRef;
	resUpdates[m_nbrEqs*iSol + evVarIDs[i]] = m_omegaTv[i]*ovOmegavRef;
      }

      // NOTE (Finding 8): pe*div(v) term for electron energy NOT implemented.
//       if (m_library->presenceElectron()) {
//         computePeDivV(element,m_srcTerm,jacob);
//       }

      CFLog(DEBUG_MAX,"ChemNEQST::computeSource() => source = " << resUpdates << "\n");

//       for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID)
//       {
//         rhs[resID] += resFactor*m_solPntJacobDets[iSol]*m_srcTerm[iEq];
//       }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTermVS::computeSourceVT(RealVector& omegaTv, CFreal& omegaRad)
{
  CFreal pdim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::P]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::P];
  CFreal Tdim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
  CFreal rhodim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO];
  m_library->getSourceTermVT(Tdim, m_tvDim, pdim, rhodim,omegaTv,omegaRad);
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTermVS::getSToStateJacobian(const CFuint iState)
{
  // Reset the Jacobian
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq] = 0.0;
  }

  // Check doComputeST (same logic as addSourceTerm)
  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint iEqSS = eqSS.getEqSS();
  SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel();
  const CFuint nbSpecies = term->getNbScalarVars(0);
  const CFuint nbEvEqs = term->getNbScalarVars(1);
  const CFuint nbEulerEq = m_dim + 2;
  const CFuint nbEqs = eqSS.getNbEqsSS();
  const vector<CFuint>& varIDs = MultiScalarVarSet<Euler2DVarSet>::EULERSET::getEqSetData()[0].getEqSetVarIDs();

  bool doComputeST = false;
  if (varIDs[0] > 0 && (iEqSS == 0 && nbEqs >= nbSpecies)) {
    doComputeST = true;
  }
  if ((varIDs[0] == 0 && (iEqSS == 0) && (nbEqs >= nbEulerEq+nbSpecies)) ||
      (varIDs[0] == 0 && (iEqSS == 1))) {
    doComputeST = true;
  }
  if (!doComputeST) return;

  static bool loggedAnalyticalPath = false;
  if (!loggedAnalyticalPath)
  {
    CFLog(INFO, "TNEQSourceTermVS::getSToStateJacobian() using analytical source Jacobian from physical-chemical library.\n");
    loggedAnalyticalPath = true;
  }

  // Compute physical data at this solution point
  State& state = *((*m_cellStates)[iState]);
  m_eulerVarSet->computePhysicalData(state, m_solPhysData);

  RealVector& refData = m_eulerVarSet->getModel()->getReferencePhysicalData();

  CFreal pdim = (m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::P] + m_eulerVarSet->getModel()->getPressInf()) *
    refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::P];
  CFreal Tdim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T] * refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
  CFreal rhodim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO] * refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO];

  const CFuint firstSpecies = term->getFirstScalarVar(0);
  for (CFuint i = 0; i < nbSpecies; ++i) {
    m_ys[i] = m_solPhysData[firstSpecies + i];
  }

  setVibTemperature(m_solPhysData, state, m_tvDim);
  m_tvDim *= refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];

  // Call PLATO with flagJac=true to get analytical Jacobian
  m_omegaTv = 0.0;
  m_omegaRad = 0.0;
  m_platoJacob = 0.0;

  if (this->m_library->getName() != "Mutation2OLD"
      && this->m_library->getName() != "MutationPanesi"
      && this->m_library->getName() != "Mutationpp") {
    this->m_library->getSource(Tdim, m_tvDim, pdim, rhodim, m_ys,
                               true, m_omega, m_omegaTv, m_omegaRad, m_platoJacob);
  }
  else {
    // Fallback: call getMassProductionTerm with flagJac=true
    RealMatrix jacob(m_nbrEqs, m_nbrEqs);
    m_library->getMassProductionTerm(Tdim, m_tvDim, pdim, rhodim, m_ys,
                                     true, m_omega, jacob);
    m_platoJacob = jacob;
  }

  // m_platoJacob(i,j) = d(prodterm[i]) / d(W_dim[j])
  // W_dim ordering: [rho_0, ..., rho_{NS-1}, (mom_0, ..., mom_{nDim-1}), T, Tv1, ...]
  //
  // For RhoivtTv state variables, the variable transformation is diagonal:
  //   state = [rho_i/refRho, ..., u/refV, v/refV, T/refT, Tv/refT]
  //   W_dim[k] = state[k] * refScale[k]
  //   where refScale = [refRho, ..., refRho, refV, ..., refV, refT, ..., refT]

  const CFreal refRho = refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO];
  const CFreal refV   = refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::V];
  const CFreal refT   = refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];

  // Build the diagonal reference scale vector using equation IDs
  // refScale[k] = dW_dim[k] / dstate[k]
  // (Finding 2 fix: use speciesVarIDs/evVarIDs instead of positional indexing)
  const vector<CFuint>& speciesVarIDs = MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].getEqSetVarIDs();
  const vector<CFuint>& evVarIDs = MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[1].getEqSetVarIDs();
  const CFuint TID = nbSpecies + m_dim;

  RealVector refScale(m_nbrEqs);
  for (CFuint i = 0; i < nbSpecies; ++i) {
    refScale[speciesVarIDs[i]] = refRho;
  }
  for (CFuint i = 0; i < m_dim; ++i) {
    refScale[nbSpecies + i] = refV;
  }
  refScale[TID] = refT;
  for (CFuint i = 0; i < nbEvEqs; ++i) {
    refScale[evVarIDs[i]] = refT;
  }

  // Source scaling factors (same as addSourceTerm)
  const CFreal ovOmegaRef = PhysicalModelStack::getActive()->getImplementor()->
    getRefLength() / (refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO] * refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::V]);
  const CFreal ovOmegavRef = PhysicalModelStack::getActive()->getImplementor()->
    getRefLength() / ((*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO] * (*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::H] * (*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::V]);

  // Chain rule: dResUpdates[iEq]/dstate[jVar] = scaling * platoJacob(platoRow, jVar) * refScale[jVar]
  // PLATO row ordering: [species_0..NS-1, (momentum), energy, ev_energy_0..]
  // COOLFluiD state ordering: uses speciesVarIDs, evVarIDs (correct for Te-enabled layouts)
  for (CFuint jVar = 0; jVar < m_nbrEqs; ++jVar)
  {
    const CFreal scale_j = refScale[jVar];

    // Species source rows (PLATO rows 0..NS-1 -> COOLFluiD speciesVarIDs)
    for (CFuint i = 0; i < nbSpecies; ++i) {
      m_stateJacobian[jVar][speciesVarIDs[i]] = ovOmegaRef * m_platoJacob(i, jVar) * scale_j;
    }

    // Total energy row (PLATO row NS+nDim -> COOLFluiD TID)
    m_stateJacobian[jVar][TID] = -ovOmegavRef * m_platoJacob(nbSpecies + m_dim, jVar) * scale_j;

    // Vibrational energy rows (PLATO rows NS+nDim+1+i -> COOLFluiD evVarIDs[i])
    for (CFuint i = 0; i < nbEvEqs; ++i) {
      m_stateJacobian[jVar][evVarIDs[i]] = ovOmegavRef * m_platoJacob(nbSpecies + m_dim + 1 + i, jVar) * scale_j;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTermVS::setVibTemperature(const RealVector& pdata,
					      const Framework::State& state,
					      RealVector& tvib)
{
  // Use equation IDs from the variable set data instead of positional
  // indexing (Finding 1 fix: correct for Te-enabled layouts where
  // vibrational temperatures may not start at species+dim+1)
  const vector<CFuint>& evVarIDs =
    Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DVarSet>::getEqSetData()[1].getEqSetVarIDs();

  for (CFuint i = 0; i < tvib.size(); ++i) {
    tvib[i] = state[evVarIDs[i]];
  }
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTermVS::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  CNEQSourceTermVS::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTermVS::setup()
{
  CFAUTOTRACE;
  CNEQSourceTermVS::setup();

  const CFuint nbVibEnergyEqs = this->m_eulerVarSet->getModel()->getNbScalarVars(1);
  m_omegaTv.resize(nbVibEnergyEqs);

  Common::SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = this->m_eulerVarSet->getModel();
  m_refData = &term->getReferencePhysicalData();

  m_platoJacob.resize(m_nbrEqs, m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTermVS::unsetup()
{
  CFAUTOTRACE;
  CNEQSourceTermVS::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    TNEQSourceTermVS::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = CNEQSourceTermVS::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
