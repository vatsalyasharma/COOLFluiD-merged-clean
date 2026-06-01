//Vatsalya: new class (dead code); inverse matrix transformer dRhoivtTv/dCons was never implemented; builds forward matrix dCons/dRhoivtTv and inverts via LU; has ZERO runtime effect because FinalizeRHSCom=Null in production

#include "NEQ/NEQ.hh"
#include "Euler2DNEQConsToRhoivtTvInRhoivtTvVS.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQConsToRhoivtTvInRhoivtTvVS,
			    VarSetMatrixTransformer,
			    NEQModule, 1>
euler2DNEQConsToRhoivtTvInRhoivtTvVSProvider("Euler2DNEQConsToRhoivtTvInRhoivtTvVS");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQConsToRhoivtTvInRhoivtTvVS::Euler2DNEQConsToRhoivtTvInRhoivtTvVS
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ys(),
  _tvDim(),
  _evDim(),
  _dhe(),
  _fwdMatrix(),
  _inverter(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQConsToRhoivtTvInRhoivtTvVS::~Euler2DNEQConsToRhoivtTvInRhoivtTvVS()
{
  delete _inverter;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQConsToRhoivtTvInRhoivtTvVS::setMatrix(const RealVector& state)
{
  cf_assert(_model.isNotNull());

  // find a way of storing this pointer (add setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  Common::SafePtr<PhysicalChemicalLibrary::ExtraData> eData = library->getExtraData();

  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint nbEqs = nbSpecies + 4; // species + u + v + T + Tv

  // Lazy initialization of inverter and forward matrix
  if (_inverter == CFNULL) {
    _inverter = MatrixInverter::create(nbEqs, false);
    _fwdMatrix.resize(nbEqs, nbEqs);
  }

  // Set the mixture density (sum of the partial densities)
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += state[ie];
  }

  // Set the species fractions
  const CFreal ovRho = 1./rho;
  _ys.resize(nbSpecies);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ys[ie] = state[ie]*ovRho;
  }

  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!!
  library->setSpeciesFractions(_ys);

  const RealVector& refData = _model->getReferencePhysicalData();
  CFreal rhoDim = rho*refData[EulerTerm::RHO];
  CFreal T = state[nbSpecies + 2];
  CFreal Tdim = T*refData[EulerTerm::T];

  // set the vibrational temperature
  const CFuint nbTv = _model->getNbScalarVars(1);
  cf_assert(nbTv == 1);
  _dhe.resize(3 + nbTv);
  _tvDim.resize(nbTv);
  _evDim.resize(nbTv);

  const CFuint startTv = nbSpecies + 3;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    _tvDim[ie] = state[startTv + ie]*refData[EulerTerm::T];
  }

  CFreal p = library->pressure(rhoDim, Tdim, &_tvDim[0]);
  CFreal pdim = p*refData[EulerTerm::P];
  const CFreal u = state[nbSpecies];
  const CFreal v = state[nbSpecies + 1];

  // vibrational temperatures
  library->setDensityEnthalpyEnergy(Tdim, _tvDim, pdim, _dhe, true);

  const CFuint uID  = nbSpecies;
  const CFuint vID  = nbSpecies+1;
  const CFuint eID  = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  const CFreal eT = eData->dEdT;
  const CFreal evTv = eData->dEvTv;
  const CFreal q = 0.5*(u*u + v*v);

  // Step 1: Build forward matrix dCons/dRhoivtTv (same logic as the inverse class)
  _fwdMatrix = 0.0;

  for (CFuint is = 0; is < nbSpecies; ++is) {
    _fwdMatrix(is,is) = 1.0;
    _fwdMatrix(uID,is) = u;
    _fwdMatrix(vID,is) = v;
    _fwdMatrix(eID,is) = eData->dRhoEdRhoi[is] + eData->dRhoEvdRhoi[is] + q;
    _fwdMatrix(evID,is) = eData->dRhoEvdRhoi[is];
  }

  _fwdMatrix(uID,uID) = rho;
  _fwdMatrix(vID,vID) = rho;

  _fwdMatrix(eID,uID)  = rho*u;
  _fwdMatrix(eID,vID)  = rho*v;
  _fwdMatrix(eID,eID)  = rho*eT;
  _fwdMatrix(eID,evID) = rho*evTv;

  _fwdMatrix(evID,evID) = rho*evTv;

  // Step 2: Invert to get dRhoivtTv/dCons
  _transMatrix = 0.0;
  _inverter->invert(_fwdMatrix, _transMatrix);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
