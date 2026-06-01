#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"
#include "MathTools/MathChecks.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/PhysicalityEuler2DVS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace {

//Vatsalya: added NaN/Inf guard helpers; original had no protection against non-finite states contaminating cell averages and poisoning all solution points
template <typename VecT>
bool hasNonFiniteEntries(const VecT& vec, const COOLFluiD::CFuint nEq)
{
  for (COOLFluiD::CFuint iEq = 0; iEq < nEq; ++iEq)
  {
    if (!cfFinite(vec[iEq]))
    {
      return true;
    }
  }
  return false;
}

template <typename VecT>
void sanitizeNonFiniteEntries(VecT& vec, const COOLFluiD::CFuint nEq, const COOLFluiD::CFreal fallback = 0.0)
{
  for (COOLFluiD::CFuint iEq = 0; iEq < nEq; ++iEq)
  {
    if (!cfFinite(vec[iEq]))
    {
      vec[iEq] = fallback;
    }
  }
}

} // namespace

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PhysicalityEuler2DVS, FluxReconstructionSolverData, FluxReconstructionNavierStokesModule>
    PhysicalityEuler2DVSFRProvider("PhysicalityEuler2DVS");

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler2DVS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("MinDensity","Minimum allowable value for density (only used for Cons).");
  options.addConfigOption< CFreal >("MinPressure","Minimum allowable value for pressure.");
  options.addConfigOption< CFreal >("MinTemperature","Minimum allowable value for temperature (only used for Puvt).");
  options.addConfigOption< bool >("CheckInternal","Boolean to tell wether to also check internal solution for physicality.");
  options.addConfigOption< bool >("LimCompleteState","Boolean to tell wether to limit complete state or single variable.");
  options.addConfigOption< bool >("ExpLim","Boolean to tell wether to use the experimental limiter.");
  options.addConfigOption< bool >("GammaClipping","Boolean to tell wether to clip gamma.");
  options.addConfigOption< std::vector<CFreal> >("MinTurbVars","Minimum K, Omega,... values");
  options.addConfigOption< CFreal >("MinGamma","Minimum clipping value for gamma.");
  options.addConfigOption< CFreal >("MaxGamma","Maximum clipping value for gamma.");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityEuler2DVS::PhysicalityEuler2DVS(const std::string& name) :
  BasePhysicality(name),
  m_minDensity(),
  m_minPressure(),
  m_minTemperature(),
  m_eulerVarSet(CFNULL),
  m_eulerVarSetMS(CFNULL),
  m_gammaMinusOne(),
  m_solPhysData(),
  m_cellAvgState(),
  m_cellAvgSolCoefs(),
  m_nbSpecies()
{
  addConfigOptionsTo(this);

  m_minDensity = 1e-2;
  setParameter( "MinDensity", &m_minDensity );

  m_minPressure = 1e-2;
  setParameter( "MinPressure", &m_minPressure );

  m_minTemperature = 1e-2;
  setParameter( "MinTemperature", &m_minTemperature );

  m_checkInternal = false;
  setParameter( "CheckInternal", &m_checkInternal );

  m_limCompleteState = true;
  setParameter( "LimCompleteState", &m_limCompleteState );

  m_expLim = true;
  setParameter( "ExpLim", &m_expLim );

  m_clipGamma = false;
  setParameter( "GammaClipping", &m_clipGamma );

  m_minTurbVars = vector<CFreal>();
  setParameter("MinTurbVars",&m_minTurbVars);

  m_minGamma = 0.01;
  setParameter( "MinGamma", &m_minGamma );

  m_maxGamma = 0.99;
  setParameter( "MaxGamma", &m_maxGamma );
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityEuler2DVS::~PhysicalityEuler2DVS()
{
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler2DVS::configure ( Config::ConfigArgs& args )
{
  BasePhysicality::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

bool PhysicalityEuler2DVS::checkPhysicality()
{
  bool physical = true;
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt" || getMethodData().getUpdateVarStr() == "BSLPuvt" || getMethodData().getUpdateVarStr() == "SSTPuvt";
  const bool Cons = getMethodData().getUpdateVarStr() == "Cons";
  const bool RhoivtLike = getMethodData().getUpdateVarStr() == "RhoivtTv" ||
                          getMethodData().getUpdateVarStr() == "Rhoivt"; //Vatsalya: original checked only "RhoivtTv" (TCNEQ 2T); CNEQ uses "Rhoivt" (1T) which bypassed all NEQ-specific gradient/composition/physicality paths
  //const bool hasArtVisc = getMethodData().hasArtificialViscosity();
  //DataHandle< CFreal > posPrev = socket_posPrev.getDataHandle();
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();
  const CFuint cellID = m_cell->getID();

  const CFuint nbTurbVars = m_cellStatesFlxPnt[0].size() - 4;

  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    if (hasNonFiniteEntries(m_cellStatesFlxPnt[iFlx], m_nbrEqs))
    {
      physical = false;
      continue;
    }

    if (Puvt)
    {
      if (m_cellStatesFlxPnt[iFlx][0] < m_minPressure || m_cellStatesFlxPnt[iFlx][3] < m_minTemperature)
      {
	physical = false;
      }

      for (CFuint iTurb = 0; iTurb < nbTurbVars; ++iTurb)
      {
        if (m_cellStatesFlxPnt[iFlx][4+iTurb] < m_minTurbVars[iTurb])
        {
	  physical = false;
        }
      }

//      if (hasArtVisc)
//      {
//	posPrev[cellID] = min(m_cellStatesFlxPnt[iFlx][0]/m_minPressure,posPrev[cellID]);
//	posPrev[cellID] = min(m_cellStatesFlxPnt[iFlx][3]/m_minTemperature,posPrev[cellID]);
//      }
    }
    else if(Cons)
    {
      CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
      CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
      CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
      CFreal rhoE = m_cellStatesFlxPnt[iFlx][3];
      if (rho < m_minDensity)
      {
        physical = false;
      }
      else
      {
        const CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
        if (!cfFinite(press) || press < m_minPressure)
        {
          physical = false;
        }
      }

//      if (hasArtVisc)
//      {
//	posPrev[cellID] = min(rho/m_minDensity,posPrev[cellID]);
//	posPrev[cellID] = min(press/m_minPressure,posPrev[cellID]);
//      }
    }
    else if(RhoivtLike)
    {
      //Vatsalya: original had hardcoded `i!=2 && i!=4` species check (5-species air assumption); replaced with generic loop over m_nbSpecies for any mixture
      for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	CFreal rho  = m_cellStatesFlxPnt[iFlx][i];
	if (rho < m_minDensity)
	{
	  physical = false;
	  break;
	}
      }
      for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	CFreal T = m_cellStatesFlxPnt[iFlx][i];
	if( T <  m_minTemperature ){
	  physical = false;
	  break;
	}
      }

    }
    //cout << "here1 Bis2"<<endl;
  }
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    output[((*m_cellStates)[iSol])->getLocalID()] = 0.0;
  }

  if (physical && m_checkInternal)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      if (hasNonFiniteEntries(*((*m_cellStates)[iSol]), m_nbrEqs))
      {
        physical = false;
        continue;
      }

	      if (Puvt)
	      {
	        if ((*((*m_cellStates)[iSol]))[0] < m_minPressure || (*((*m_cellStates)[iSol]))[3] < m_minTemperature)
	        {
	  physical = false;
        }

        for (CFuint iTurb = 0; iTurb < nbTurbVars; ++iTurb)
        {
          if ((*((*m_cellStates)[iSol]))[4+iTurb] < m_minTurbVars[iTurb])
          {
	    physical = false;
          }
        }
      }
      else if(Cons)
      {
        CFreal rho  = (*((*m_cellStates)[iSol]))[0];

	if (rho < m_minDensity)
        {
	  physical = false;
        }
        else
	{
	  CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
          CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
	      CFreal rhoE = (*((*m_cellStates)[iSol]))[3];
		  CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);

		  if (!cfFinite(press) || press < m_minPressure)
		  {
		    physical = false;
		  }
	}
      }
      else if (RhoivtLike){
	//Vatsalya: original had hardcoded `i!=2 && i!=4` species check (5-species air assumption); replaced with generic loop over m_nbSpecies for any mixture
	for (CFuint i = 0 ; i< m_nbSpecies ; ++i){
	  if((*((*m_cellStates)[iSol]))[i] < m_minDensity){
	    physical = false;
	    break;
	  }
	}
	for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	  if((*((*m_cellStates)[iSol]))[i] < m_minTemperature){
	    physical = false;
	    break;
	  }
	}
      }
    }
  }

  // clip gamma if needed (default false, only activate in GReKO or gamma-alpha)
  if (Puvt && m_clipGamma)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      if ((*((*m_cellStates)[iSol]))[6] < m_minGamma)
      {
	(*((*m_cellStates)[iSol]))[6] = m_minGamma;
      }
      else if ((*((*m_cellStates)[iSol]))[6] > m_maxGamma)
      {
        (*((*m_cellStates)[iSol]))[6] = m_maxGamma;
      }
    }
  }

  return physical;

}
//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler2DVS::enforcePhysicality()
{
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt" || getMethodData().getUpdateVarStr() == "BSLPuvt" || getMethodData().getUpdateVarStr() == "SSTPuvt";
  const bool Cons = getMethodData().getUpdateVarStr() == "Cons";
  const bool RhoivtLike = getMethodData().getUpdateVarStr() == "RhoivtTv" ||
                          getMethodData().getUpdateVarStr() == "Rhoivt";
  bool needsLim = false;
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();

  const CFuint nbTurbVars = m_cellStatesFlxPnt[0].size() - 4;

  // compute average state
  m_cellAvgState = 0.;
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellAvgState += (*m_cellAvgSolCoefs)[iSol] * (*(*m_cellStates)[iSol]);
  }

  if (hasNonFiniteEntries(m_cellAvgState, m_nbrEqs))
  {
    needsLim = true;
  }

  // check if average state is physical
  if (Puvt)
  {
    if (m_cellAvgState[0] < m_minPressure || m_cellAvgState[3] < m_minTemperature) needsLim = true;

    for (CFuint iTurb = 0; iTurb < nbTurbVars; ++iTurb)
    {
      if (m_cellAvgState[4+iTurb] < m_minTurbVars[iTurb])
      {
        needsLim = true;
      }
    }
  }
  else if (Cons)
  {
    CFreal rho  = m_cellAvgState[0];

    if (rho < m_minDensity)
    {
      needsLim = true;
    }
    else
    {
      CFreal rhoU = m_cellAvgState[1];
      CFreal rhoV = m_cellAvgState[2];
      CFreal rhoE = m_cellAvgState[3];
      CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);

      if (!cfFinite(press) || press < m_minPressure) needsLim = true;
    }
  }
  else if (RhoivtLike){
    for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
      if(m_cellAvgState[i] < m_minDensity){
	needsLim = true;
      }
    }
    for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
      if(m_cellAvgState[i] < m_minTemperature){
	needsLim = true;
      }
    }
  }

  // if average state is unphysical, modify the unphysical variable
  if (needsLim)
  {
    m_nbAvLimits += 1;

    // Prevent NaN/Inf average entries from contaminating all solution points.
    sanitizeNonFiniteEntries(m_cellAvgState, m_nbrEqs);

    // subtract the average solution from the state
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      output[((*m_cellStates)[iSol])->getLocalID()] = -10000.0;
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
	(*((*m_cellStates)[iSol]))[iEq] -= m_cellAvgState[iEq];
      }
    }

    if (Puvt)
    {
      if (m_cellAvgState[0] < m_minPressure)
      {
        m_cellAvgState[0] = 1.1*m_minPressure;
      }
      if (m_cellAvgState[3] < m_minTemperature)
      {
	m_cellAvgState[3] = 1.1*m_minTemperature;
      }

      for (CFuint iTurb = 0; iTurb < nbTurbVars; ++iTurb)
      {
        if (m_cellAvgState[4+iTurb] < m_minTurbVars[iTurb])
        {
          m_cellAvgState[4+iTurb] = 1.1*m_minTurbVars[iTurb];
        }
      }
    }
    else if (Cons)
    {
      if (m_cellAvgState[0] < m_minDensity)
      {
        m_cellAvgState[0] = 1.1*m_minDensity;
      }

      CFreal rho = m_cellAvgState[0];
      CFreal rhoU = m_cellAvgState[1];
      CFreal rhoV = m_cellAvgState[2];
      CFreal rhoE = m_cellAvgState[3];
      CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);

      if (!cfFinite(press) || press < m_minPressure)
      {
	m_cellAvgState[3] = 1.1*m_minPressure/m_gammaMinusOne + 0.5*(rhoU*rhoU+rhoV*rhoV)/rho;
      }
    }
    else if (RhoivtLike){
      for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	if(m_cellAvgState[i] < m_minDensity){
	  m_cellAvgState[i] = 1.1*m_minDensity;
	}
      }
      for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	if(m_cellAvgState[i] < m_minTemperature){
	  m_cellAvgState[i] = 1.1*m_minTemperature;
	}
      }
    }

    // compute the new states
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
	(*((*m_cellStates)[iSol]))[iEq] += m_cellAvgState[iEq];
      }
    }

    // compute the new states in the flux points
    computeFlxPntStates(m_cellStatesFlxPnt);
  }

  // flags telling which state needs to be limited
  vector<bool> needsLimFlags(m_nbrEqs);

  if (!m_expLim)
  {

  // check each flux point for unphysical states
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {


    // reset limit flags
    needsLim = false;
    for (CFuint iFlag = 0; iFlag < m_nbrEqs; ++iFlag)
    {
      needsLimFlags[iFlag] = false;
    }

    if (Puvt)
    {
      if (m_cellStatesFlxPnt[iFlx][0] < m_minPressure || m_cellStatesFlxPnt[iFlx][3] < m_minTemperature)
      {
	needsLim = true;
	needsLimFlags[0] = m_cellStatesFlxPnt[iFlx][0] < m_minPressure;
	needsLimFlags[3] = m_cellStatesFlxPnt[iFlx][3] < m_minTemperature;
      }
    }
    else if (Cons)
    {
      CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
      CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
      CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
      CFreal rhoE = m_cellStatesFlxPnt[iFlx][3];
      CFreal press = 0.0;
      if (rho >= m_minDensity)
      {
        press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
      }

      if (rho < m_minDensity || !cfFinite(press) || press < m_minPressure)
      {
	needsLim = true;
	needsLimFlags[0] = rho < m_minDensity;
	needsLimFlags[3] = (!cfFinite(press) || press < m_minPressure);
      }
    }
    else if (RhoivtLike){
      for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	if (m_cellStatesFlxPnt[iFlx][i] < m_minDensity){
	  needsLim = true;
	  needsLimFlags[i] = true;
	}
      }
      for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	if(m_cellStatesFlxPnt[iFlx][i] < m_minTemperature){
	  needsLim = true;
	  needsLimFlags[i] = true;
	}
      }
    }

    // if needed, limit states
    if (needsLim)
    {
      // limiting factor
      CFreal phi = 1.0;

      for (CFuint iScale = 0; iScale < 10; ++iScale)
      {
        phi /= 2.0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	  {
	    if( needsLimFlags[iEq] || m_limCompleteState )
	    {
	      (*((*m_cellStates)[iSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*((*((*m_cellStates)[iSol]))[iEq]);
	    }
	  }
        }

        // recompute the states in the flux points
        computeFlxPntStates(m_cellStatesFlxPnt);

	// reset needsLim
	needsLim = false;

	// reset limit flags
        for (CFuint iFlag = 0; iFlag < m_nbrEqs; ++iFlag)
        {
          needsLimFlags[iFlag] = false;
        }

	// check if the state is now physical
		if (hasNonFiniteEntries(m_cellStatesFlxPnt[iFlx], m_nbrEqs))
                {
                  needsLim = true;
                  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                  {
                    needsLimFlags[iEq] = !cfFinite(m_cellStatesFlxPnt[iFlx][iEq]);
                  }
                }
		else if (Puvt)
	        {
	          if (m_cellStatesFlxPnt[iFlx][0] < m_minPressure || m_cellStatesFlxPnt[iFlx][3] < m_minTemperature)
	          {
	    needsLim = true;
	    needsLimFlags[0] = m_cellStatesFlxPnt[iFlx][0] < m_minPressure;
	    needsLimFlags[3] = m_cellStatesFlxPnt[iFlx][3] < m_minTemperature;
          }
        }
	        else if (Cons)
	        {
	          CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
	          CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
	          CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
	          CFreal rhoE = m_cellStatesFlxPnt[iFlx][3];
	          CFreal press = 0.0;
                  if (rho >= m_minDensity)
                  {
                    press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
                  }

	          if (rho < m_minDensity || !cfFinite(press) || press < m_minPressure)
	          {
		    needsLim = true;
		    needsLimFlags[0] = rho < m_minDensity;
		    needsLimFlags[3] = (!cfFinite(press) || press < m_minPressure);
	          }
	        }
	else if (RhoivtLike){
	  for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	    if (m_cellStatesFlxPnt[iFlx][i] < m_minDensity){
	      needsLim = true;
	      needsLimFlags[i] = true;
	    }
	  }
	  for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	    if(m_cellStatesFlxPnt[iFlx][i] < m_minTemperature){
	      needsLim = true;
	      needsLimFlags[i] = true;
	    }
	  }
	}
	// break if the states are physical
	if(!needsLim)
	{
	  break;
	}
      }

      // if after limiting the states are still non-physical, set them to the average states
      if (needsLim)
	{
	  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	    {
	      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
		{
		  if ( needsLimFlags[iEq] || m_limCompleteState)
		  {
		    (*((*m_cellStates)[iSol]))[iEq] = m_cellAvgState[iEq];
		    //cout << " give the av sol  "<< m_cellAvgState[iEq] << endl;
		  }
		}
	    }
	    break;
	}
    }
  }
  }
  else
  {
    bool nonFiniteFlx = false;
    for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
    {
      if (hasNonFiniteEntries(m_cellStatesFlxPnt[iFlx], m_nbrEqs))
      {
        nonFiniteFlx = true;
        break;
      }
    }
    if (nonFiniteFlx)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*((*m_cellStates)[iSol]))[iEq] = m_cellAvgState[iEq];
        }
      }
      computeFlxPntStates(m_cellStatesFlxPnt);
    }

    if (Cons)
    {
      CFreal rhoAv = m_cellAvgState[0];
      if (rhoAv < m_minDensity)
      {
        rhoAv = m_minDensity;
      }
      CFreal rhoUAv = m_cellAvgState[1];
      CFreal rhoVAv = m_cellAvgState[2];
      CFreal rhoEAv = m_cellAvgState[3];
      CFreal pressAv = m_gammaMinusOne*(rhoEAv - 0.5*(rhoUAv*rhoUAv+rhoVAv*rhoVAv)/rhoAv);
      if (!cfFinite(pressAv))
      {
        pressAv = m_minPressure;
      }

      CFreal epsilon = min(rhoAv,pressAv);
      epsilon = min(m_minDensity,epsilon);
      //CFreal epsilonP = 0.8*pressAv;
      CFreal epsilonP = min(m_minPressure,epsilon);


      CFreal rhoMin = 1.0e13;

      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
      {
        rhoMin = min(rhoMin,m_cellStatesFlxPnt[iFlx][0]);
      }

      CFreal coeff = min((rhoAv-epsilon)/(rhoAv-rhoMin),1.0);

      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[0] = (1.0-coeff)*m_cellAvgState[0] + coeff*((*((*m_cellStates)[iSol]))[0]);
        }
      }

      // recompute the states in the flux points
      computeFlxPntStates(m_cellStatesFlxPnt);

      CFreal t = 1.0;

	      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
	      {
	        CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
	        CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
	        CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
	        CFreal rhoE = m_cellStatesFlxPnt[iFlx][3];
                if (rho <= m_minDensity)
                {
                  t = 0.0;
                  continue;
                }
	        CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);

		//CFLog(INFO, "p " << iFlx << " : "<< press << "\n");

	        if (!cfFinite(press))
                {
                  t = 0.0;
                  continue;
                }
                if (press < epsilon)
	        {
		  CFreal A = rho*rhoE+rhoAv*rhoEAv-rhoAv*rhoE-rho*rhoEAv-0.5*(rhoUAv*rhoUAv+rhoU*rhoU-2.0*rhoU*rhoUAv+rhoVAv*rhoVAv+rhoV*rhoV-2.0*rhoV*rhoVAv);
		  CFreal B = rhoAv*rhoE+rho*rhoEAv-2.0*rhoAv*rhoEAv-epsilonP/m_gammaMinusOne*(rho-rhoAv)-0.5*(2.0*rhoU*rhoUAv-2.0*rhoUAv*rhoUAv+2.0*rhoV*rhoVAv-2.0*rhoVAv*rhoVAv);
		  CFreal C = rhoAv*rhoEAv-0.5*(rhoUAv*rhoUAv+rhoVAv*rhoVAv)-epsilonP/m_gammaMinusOne*rhoAv;
		  CFreal D = B*B-4.0*A*C;
                  if (D < 0.0 || !cfFinite(D))
                  {
                    t = 0.0;
                    continue;
                  }
		  CFreal sol1 = (-B+sqrt(D))/(2.0*A);
	  if (sol1 < 0.0 || sol1 > 1.0)
	  {
	    sol1 = (-B-sqrt(D))/(2.0*A);
	  }
	  cf_assert(sol1>-1.0e-13 && sol1<1.000001);

	  t = min(t,sol1);
        }
      }
      ////////
      if (t < 1.0)
      {
	//CFLog(INFO, "t: " << t << "\n");
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  CFLog(VERBOSE, "state " << iSol << " : "<< (*((*m_cellStates)[iSol])) << ", coord: " << (*((*m_cellStates)[iSol])).getCoordinates() << "\n");
	  CFreal p = m_gammaMinusOne*((*((*m_cellStates)[iSol]))[3] - 0.5*(pow((*((*m_cellStates)[iSol]))[1],2)+pow((*((*m_cellStates)[iSol]))[2],2))/(*((*m_cellStates)[iSol]))[0]);
	  CFLog(VERBOSE, "p " << iSol << " : "<< p << "\n");
          CFLog(VERBOSE, "T " << iSol << " : "<< p/((*((*m_cellStates)[iSol]))[0]) << " (p/rho)\n"); //Vatsalya: original used hardcoded R=287.046 (wrong for non-air); replaced with p/rho ratio
	  CFLog(VERBOSE, "u " << iSol << " : "<< (*((*m_cellStates)[iSol]))[1]/(*((*m_cellStates)[iSol]))[0] << ", v " << iSol << " : "<< (*((*m_cellStates)[iSol]))[2]/(*((*m_cellStates)[iSol]))[0] << "\n");
	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
          {
            (*((*m_cellStates)[iSol]))[iEq] = (1.0-t)*m_cellAvgState[iEq] + t*((*((*m_cellStates)[iSol]))[iEq]);
          }
          CFLog(VERBOSE, "state2 " << iSol << " : "<< (*((*m_cellStates)[iSol])) << "\n");
	  p = m_gammaMinusOne*((*((*m_cellStates)[iSol]))[3] - 0.5*(pow((*((*m_cellStates)[iSol]))[1],2)+pow((*((*m_cellStates)[iSol]))[2],2))/(*((*m_cellStates)[iSol]))[0]);
	  CFLog(VERBOSE, "p2 " << iSol << " : "<< p << "\n");
          CFLog(VERBOSE, "T2 " << iSol << " : "<< p/((*((*m_cellStates)[iSol]))[0]) << " (p/rho)\n"); //Vatsalya: original used hardcoded R=287.046 (wrong for non-air); replaced with p/rho ratio
	  CFLog(VERBOSE, "u2 " << iSol << " : "<< (*((*m_cellStates)[iSol]))[1]/(*((*m_cellStates)[iSol]))[0] << ", v2 " << iSol << " : "<< (*((*m_cellStates)[iSol]))[2]/(*((*m_cellStates)[iSol]))[0] << "\n");
        }

      }
    }
    else if (Puvt)
    {
      CFreal pAv = m_cellAvgState[0];
      CFreal TAv = m_cellAvgState[3];

      CFreal pMin = 1.0e13;

      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
      {
        pMin = min(pMin,m_cellStatesFlxPnt[iFlx][0]);
      }

      CFreal coeff = min((pAv-m_minPressure)/(pAv-pMin),1.0);

      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 1.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[0] = (1.0-coeff)*m_cellAvgState[0] + coeff*((*((*m_cellStates)[iSol]))[0]);
        }
      }

      CFreal TMin = 1.0e13;

      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
      {
        TMin = min(TMin,m_cellStatesFlxPnt[iFlx][3]);
      }

      coeff = min((TAv-m_minTemperature)/(TAv-TMin),1.0);

      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[3] = (1.0-coeff)*m_cellAvgState[3] + coeff*((*((*m_cellStates)[iSol]))[3]);
        }
      }

      for (CFuint iTurb = 0; iTurb < nbTurbVars; ++iTurb)
      {
        CFreal turbMin = 1.0e13;

        CFreal turbAv = m_cellAvgState[iTurb+4];

        for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
        {
          turbMin = min(turbMin,m_cellStatesFlxPnt[iFlx][iTurb+4]);
        }

        coeff = min((turbAv-m_minTurbVars[iTurb])/(turbAv-turbMin),1.0);

        if (coeff < 1.0)
        {
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
	    output[((*m_cellStates)[iSol])->getLocalID()] += pow(10.0,iTurb+2);

            (*((*m_cellStates)[iSol]))[iTurb+4] = (1.0-coeff)*turbAv + coeff*((*((*m_cellStates)[iSol]))[iTurb+4]);
          }
        }
      }
    }
    else if (RhoivtLike)
    {
      //Vatsalya: rewrote ~230 lines of buggy conservative limiter to ~70 lines of coupled Zhang-Shu; original had hardcoded R=287 (wrong for TCNEQ), hardcoded i!=2&&i!=4 (broke for 5+ species), used wrong flux point for all sol points, temperature limiting was commented out; for Rhoivt rho>0 and T>0 guarantees p>0 so no conservative pressure limiting needed
      // Coupled Zhang-Shu positivity limiter for multi-species TCNEQ
      // Two-step approach: (1) density limiting, (2) temperature limiting
      // Each step uses a SINGLE theta across ALL equations to preserve
      // polynomial shape uniformly (coupled limiting).
      // For Rhoivt/RhoivtTv, p = rho * R_mix * T, so rho>0 and T>0 => p>0.
      // No conservative pressure limiting needed.

      // Step 1: Coupled density limiting
      // Find minimum theta across all species and all flux points
      CFreal thetaRho = 1.0;
      for (CFuint i = 0; i < m_nbSpecies; ++i)
      {
        for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
        {
          if (m_cellStatesFlxPnt[iFlx][i] < m_minDensity && m_cellAvgState[i] > m_minDensity)
          {
            const CFreal t = (m_cellAvgState[i] - m_minDensity) /
                             (m_cellAvgState[i] - m_cellStatesFlxPnt[iFlx][i]);
            thetaRho = min(thetaRho, t);
          }
        }
      }

      if (thetaRho < 1.0)
      {
        // Blend ALL equations toward cell average with single theta
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
          for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
          {
            (*((*m_cellStates)[iSol]))[iEq] = (1.0-thetaRho)*m_cellAvgState[iEq] +
                thetaRho*((*((*m_cellStates)[iSol]))[iEq]);
          }
        }
        // Recompute flux point states after density limiting
        computeFlxPntStates(m_cellStatesFlxPnt);
      }

      // Step 2: Coupled temperature limiting
      // Find minimum theta across all temperatures and all flux points
      CFreal thetaT = 1.0;
      for (CFuint i = m_nbSpecies + nbDims; i < m_nbrEqs; ++i)
      {
        for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
        {
          if (m_cellStatesFlxPnt[iFlx][i] < m_minTemperature && m_cellAvgState[i] > m_minTemperature)
          {
            const CFreal t = (m_cellAvgState[i] - m_minTemperature) /
                             (m_cellAvgState[i] - m_cellStatesFlxPnt[iFlx][i]);
            thetaT = min(thetaT, t);
          }
        }
      }

      if (thetaT < 1.0)
      {
        // Blend ALL equations toward cell average with single theta
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
          for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
          {
            (*((*m_cellStates)[iSol]))[iEq] = (1.0-thetaT)*m_cellAvgState[iEq] +
                thetaT*((*((*m_cellStates)[iSol]))[iEq]);
          }
        }
      }
    }
  }


  if (m_checkInternal)
  {
    if (!m_expLim)
    {
    // chech if the solution point states are physical
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      // reset the limit flag
      needsLim = false;

      // reset limit flags
      for (CFuint iFlag = 0; iFlag < m_nbrEqs; ++iFlag)
      {
        needsLimFlags[iFlag] = false;
      }

      if (Puvt)
      {
        if ((*((*m_cellStates)[iSol]))[0] < m_minPressure || (*((*m_cellStates)[iSol]))[3] < m_minTemperature)
        {
	  needsLim = true;
	  needsLimFlags[0] = (*((*m_cellStates)[iSol]))[0] < m_minPressure;
	  needsLimFlags[3] = (*((*m_cellStates)[iSol]))[3] < m_minTemperature;
        }
      }
	      else if (Cons)
	      {
	        CFreal rho  = (*((*m_cellStates)[iSol]))[0];
	        CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
	        CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
	        CFreal rhoE = (*((*m_cellStates)[iSol]))[3];
	        CFreal press = 0.0;
                if (rho >= m_minDensity)
                {
                  press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
                }

	        if (rho < m_minDensity || !cfFinite(press) || press < m_minPressure)
	        {
		  needsLim = true;
		  needsLimFlags[0] = rho < m_minDensity;
		  needsLimFlags[3] = (!cfFinite(press) || press < m_minPressure);
	        }
	      }
      else if (RhoivtLike){
	for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	  if ((*((*m_cellStates)[iSol]))[i] < m_minDensity){
	    needsLim = true;
	    needsLimFlags[i] = true;
	  }
	}
	for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	  if((*((*m_cellStates)[iSol]))[i] <  m_minTemperature ){
	    needsLim = true;
	    needsLimFlags[i] = true;
	  }
	}
      }

      // limit the states if needed
      if (needsLim)
      {
        //CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
	// limiting factor
	CFreal phi = 1.0;

        for (CFuint iScale = 0; iScale < 10; ++iScale)
        {
          phi /= 2.0;
          for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
          {
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
	      if( needsLimFlags[iEq] || m_limCompleteState )
	      {
	        (*((*m_cellStates)[jSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*(*((*m_cellStates)[jSol]))[iEq];
	      }
	    }
          }

          // reset limit flag
          needsLim = false;
	  for (CFuint iFlag = 0; iFlag < m_nbrEqs; ++iFlag)
          {
            needsLimFlags[iFlag] = false;
          }

	  // check if the states are now physical
		  if (hasNonFiniteEntries(*((*m_cellStates)[iSol]), m_nbrEqs))
                  {
                    needsLim = true;
                    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                    {
                      needsLimFlags[iEq] = !cfFinite((*((*m_cellStates)[iSol]))[iEq]);
                    }
                  }
		  else if (Puvt)
	          {
	            if ((*((*m_cellStates)[iSol]))[0] < m_minPressure || (*((*m_cellStates)[iSol]))[3] < m_minTemperature)
	            {
	      needsLim = true;
	      needsLimFlags[0] = (*((*m_cellStates)[iSol]))[0] < m_minPressure;
	      needsLimFlags[3] = (*((*m_cellStates)[iSol]))[3] < m_minTemperature;
            }
          }
	          else if (Cons)
	          {
	            CFreal rho  = (*((*m_cellStates)[iSol]))[0];
	            CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
	            CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
	            CFreal rhoE = (*((*m_cellStates)[iSol]))[3];
	            CFreal press = 0.0;
                    if (rho >= m_minDensity)
                    {
                      press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
                    }

	            if (rho < m_minDensity || !cfFinite(press) || press < m_minPressure)
	            {
		      needsLim = true;
		      needsLimFlags[0] = rho < m_minDensity;
		      needsLimFlags[3] = (!cfFinite(press) || press < m_minPressure);
	            }
	          }
	  else if (RhoivtLike){
	    for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	      if ((*((*m_cellStates)[iSol]))[i] < m_minDensity){
		needsLim = true;
		needsLimFlags[i] = true;
	      }
	    }
	    for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	      if((*((*m_cellStates)[iSol]))[i] <  m_minTemperature ){
		needsLim = true;
		needsLimFlags[i] = true;
	      }
	    }
	  }

	  // break if the states are physical
	  if(!needsLim)
	  {
	    break;
	  }
        }

        // if still not physical after limiting, set the states to the average states
        if (needsLim)
        {
	  for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
          {
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
	      if ( needsLimFlags[iEq] || m_limCompleteState )
	      {
	        (*((*m_cellStates)[jSol]))[iEq] = m_cellAvgState[iEq];
	      }
	    }
          }
          break;
        }
      }
    }
  }
  else
  {
    bool nonFiniteSol = false;
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      if (hasNonFiniteEntries(*((*m_cellStates)[iSol]), m_nbrEqs))
      {
        nonFiniteSol = true;
        break;
      }
    }
    if (nonFiniteSol)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*((*m_cellStates)[iSol]))[iEq] = m_cellAvgState[iEq];
        }
      }
    }

    if (Cons)
    {
      CFreal rhoAv = m_cellAvgState[0];
      if (rhoAv < m_minDensity)
      {
        rhoAv = m_minDensity;
      }
      CFreal rhoUAv = m_cellAvgState[1];
      CFreal rhoVAv = m_cellAvgState[2];
      CFreal rhoEAv = m_cellAvgState[3];
      CFreal pressAv = m_gammaMinusOne*(rhoEAv - 0.5*(rhoUAv*rhoUAv+rhoVAv*rhoVAv)/rhoAv);
      if (!cfFinite(pressAv))
      {
        pressAv = m_minPressure;
      }

      CFreal epsilon = min(rhoAv,pressAv);
      epsilon = min(m_minDensity,epsilon);
      //CFreal epsilonP = 0.8*pressAv;
      CFreal epsilonP = min(m_minPressure,epsilon);


      CFreal rhoMin = 1.0e13;

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        rhoMin = min(rhoMin,(*((*m_cellStates)[iSol]))[0]);
      }

      CFreal coeff = min((rhoAv-epsilon)/(rhoAv-rhoMin),1.0);

      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[0] = (1.0-coeff)*m_cellAvgState[0] + coeff*((*((*m_cellStates)[iSol]))[0]);
        }
      }

      CFreal t = 1.0;

	      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	      {
	        CFreal rho  = (*((*m_cellStates)[iSol]))[0];
	        CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
	        CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
	        CFreal rhoE = (*((*m_cellStates)[iSol]))[3];
                if (rho <= m_minDensity)
                {
                  t = 0.0;
                  continue;
                }
	        CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);

		//CFLog(INFO, "p " << iFlx << " : "<< press << "\n");

	        if (!cfFinite(press))
                {
                  t = 0.0;
                  continue;
                }
                if (press < epsilon)
	        {
		  CFreal A = rho*rhoE+rhoAv*rhoEAv-rhoAv*rhoE-rho*rhoEAv-0.5*(rhoUAv*rhoUAv+rhoU*rhoU-2.0*rhoU*rhoUAv+rhoVAv*rhoVAv+rhoV*rhoV-2.0*rhoV*rhoVAv);
		  CFreal B = rhoAv*rhoE+rho*rhoEAv-2.0*rhoAv*rhoEAv-epsilonP/m_gammaMinusOne*(rho-rhoAv)-0.5*(2.0*rhoU*rhoUAv-2.0*rhoUAv*rhoUAv+2.0*rhoV*rhoVAv-2.0*rhoVAv*rhoVAv);
		  CFreal C = rhoAv*rhoEAv-0.5*(rhoUAv*rhoUAv+rhoVAv*rhoVAv)-epsilonP/m_gammaMinusOne*rhoAv;
		  CFreal D = B*B-4.0*A*C;
                  if (D < 0.0 || !cfFinite(D))
                  {
                    t = 0.0;
                    continue;
                  }
		  CFreal sol1 = (-B+sqrt(D))/(2.0*A);
	  if (sol1 < 0.0 || sol1 > 1.0)
	  {
	    sol1 = (-B-sqrt(D))/(2.0*A);
	  }
	  cf_assert(sol1>-1.0e-13 && sol1<1.000001);

	  t = min(t,sol1);
        }
      }
      ////////
      if (t < 1.0)
      {
	//CFLog(INFO, "t: " << t << "\n");
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  CFLog(VERBOSE, "state " << iSol << " : "<< (*((*m_cellStates)[iSol])) << ", coord: " << (*((*m_cellStates)[iSol])).getCoordinates() << "\n");
	  CFreal p = m_gammaMinusOne*((*((*m_cellStates)[iSol]))[3] - 0.5*(pow((*((*m_cellStates)[iSol]))[1],2)+pow((*((*m_cellStates)[iSol]))[2],2))/(*((*m_cellStates)[iSol]))[0]);
	  CFLog(VERBOSE, "p " << iSol << " : "<< p << "\n");
          CFLog(VERBOSE, "T " << iSol << " : "<< p/((*((*m_cellStates)[iSol]))[0]) << " (p/rho)\n"); //Vatsalya: original used hardcoded R=287.046 (wrong for non-air); replaced with p/rho ratio
	  CFLog(VERBOSE, "u " << iSol << " : "<< (*((*m_cellStates)[iSol]))[1]/(*((*m_cellStates)[iSol]))[0] << ", v " << iSol << " : "<< (*((*m_cellStates)[iSol]))[2]/(*((*m_cellStates)[iSol]))[0] << "\n");
	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
          {
            (*((*m_cellStates)[iSol]))[iEq] = (1.0-t)*m_cellAvgState[iEq] + t*((*((*m_cellStates)[iSol]))[iEq]);
          }
          CFLog(VERBOSE, "state2 " << iSol << " : "<< (*((*m_cellStates)[iSol])) << "\n");
	  p = m_gammaMinusOne*((*((*m_cellStates)[iSol]))[3] - 0.5*(pow((*((*m_cellStates)[iSol]))[1],2)+pow((*((*m_cellStates)[iSol]))[2],2))/(*((*m_cellStates)[iSol]))[0]);
	  CFLog(VERBOSE, "p2 " << iSol << " : "<< p << "\n");
          CFLog(VERBOSE, "T2 " << iSol << " : "<< p/((*((*m_cellStates)[iSol]))[0]) << " (p/rho)\n"); //Vatsalya: original used hardcoded R=287.046 (wrong for non-air); replaced with p/rho ratio
	  CFLog(VERBOSE, "u2 " << iSol << " : "<< (*((*m_cellStates)[iSol]))[1]/(*((*m_cellStates)[iSol]))[0] << ", v2 " << iSol << " : "<< (*((*m_cellStates)[iSol]))[2]/(*((*m_cellStates)[iSol]))[0] << "\n");
        }

      }
    }
    else if (Puvt)
    {
      CFreal pAv = m_cellAvgState[0];
      CFreal TAv = m_cellAvgState[3];

      CFreal pMin = 1.0e13;

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        pMin = min(pMin,(*((*m_cellStates)[iSol]))[0]);
      }

      CFreal coeff = min((pAv-m_minPressure)/(pAv-pMin),1.0);

      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 1.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[0] = (1.0-coeff)*m_cellAvgState[0] + coeff*((*((*m_cellStates)[iSol]))[0]);
        }
      }

      CFreal TMin = 1.0e13;

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        TMin = min(TMin,(*((*m_cellStates)[iSol]))[3]);
      }

      coeff = min((TAv-m_minTemperature)/(TAv-TMin),1.0);

      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[3] = (1.0-coeff)*m_cellAvgState[3] + coeff*((*((*m_cellStates)[iSol]))[3]);
        }
      }

      for (CFuint iTurb = 0; iTurb < nbTurbVars; ++iTurb)
      {
        CFreal turbMin = 1.0e13;

        CFreal turbAv = m_cellAvgState[iTurb+4];

        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          turbMin = min(turbMin,(*((*m_cellStates)[iSol]))[iTurb+4]);
        }

        coeff = min((turbAv-m_minTurbVars[iTurb])/(turbAv-turbMin),1.0);

        if (coeff < 1.0)
        {
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
	    output[((*m_cellStates)[iSol])->getLocalID()] += pow(10.0,iTurb+2);

            (*((*m_cellStates)[iSol]))[iTurb+4] = (1.0-coeff)*turbAv + coeff*((*((*m_cellStates)[iSol]))[iTurb+4]);
          }
        }
      }
    }
    else if (RhoivtLike)
    {
      CFreal rhoAvMin = 1.0e13;
      for (CFuint i = 0 ; i < m_nbSpecies ; ++i)
      {
        rhoAvMin = min(rhoAvMin,m_cellAvgState[i]);
      }

      CFreal TAvMin = 1.0e13;
      for(CFuint i = m_nbSpecies+nbDims ; i < m_nbrEqs; ++i)
      {
        TAvMin = min(TAvMin,m_cellAvgState[i]);
      }

      CFreal epsilon = min(rhoAvMin,TAvMin);
      epsilon = min(m_minDensity,epsilon);
      CFreal epsilonT = min(TAvMin,m_minTemperature);

      for (CFuint i = 0 ; i < m_nbSpecies ; ++i)
      {
        CFreal rhoMin = 1.0e13;

        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          rhoMin = min(rhoMin,(*((*m_cellStates)[iSol]))[i]);
        }

        CFreal coeff = min((m_cellAvgState[i]-epsilon)/(m_cellAvgState[i]-rhoMin),1.0);

        if (coeff < 1.0)
        {
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            (*((*m_cellStates)[iSol]))[i] = (1.0-coeff)*m_cellAvgState[i] + coeff*((*((*m_cellStates)[iSol]))[i]);
	    output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
          }
          CFLog(VERBOSE, "Lim rho " << i << "\n");
        }
      }



      for (CFuint i = m_nbSpecies+nbDims   ; i < m_nbrEqs ; ++i)
      {
        CFreal TMin = 1.0e13;

        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          TMin = min(TMin,(*((*m_cellStates)[iSol]))[i]);
        }

        CFreal coeff = min((m_cellAvgState[i]-epsilonT)/(m_cellAvgState[i]-TMin),1.0);

        if (coeff < 1.0)
        {
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            (*((*m_cellStates)[iSol]))[i] = (1.0-coeff)*m_cellAvgState[i] + coeff*((*((*m_cellStates)[iSol]))[i]);
	    output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
          }
          CFLog(VERBOSE, "Lim T " << i << "\n");
        }
      }
    }
  }
  }

//   // only needed to plot the physicality check!!!!!
//   for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//   {
//     (*((*m_cellStates)[iSol]))[0] = 10.0;
//   }
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler2DVS::setup()
{
  CFAUTOTRACE;

  // setup parent class
  BasePhysicality::setup();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const bool RhoivtLike = getMethodData().getUpdateVarStr() == "RhoivtTv" ||
                          getMethodData().getUpdateVarStr() == "Rhoivt";

  m_cellAvgSolCoefs = frLocalData[0]->getCellAvgSolCoefs();

  // get Euler 2D varset
if(!RhoivtLike){
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in PhysicalityEuler2DVSFluxReconstruction!");
  }

  // get gamma-1
  m_gammaMinusOne = m_eulerVarSet->getModel()->getGamma()-1.0;
  //m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);
}
else{
  m_eulerVarSetMS = PhysicalModelStack::getActive()-> getImplementor()->getConvectiveTerm().d_castTo< MultiScalarTerm< EulerTerm > >();
  if (m_eulerVarSetMS.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSetMS in PhysicalityEuler2DVSFluxReconstruction!");
  }
  m_nbSpecies = m_eulerVarSetMS->getNbScalarVars(0);
  // Rhoivt/RhoivtTv positivity checks use rho_i and temperature, not a calorically-perfect gamma.
  // Keep this consistent with the active model instead of hard-coding 1.4-1.
  m_gammaMinusOne = m_eulerVarSetMS->getGamma() - 1.0;
}
  m_cellAvgState.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler2DVS::unsetup()
{
  CFAUTOTRACE;

  // unsetup parent class
  BasePhysicality::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
