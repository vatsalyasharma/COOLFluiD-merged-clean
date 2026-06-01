#include "ArcJet/ArcJet.hh"
#include "ArcJet/ArcJetLTEPhysicalModel.hh"
#include "NavierStokes/NSTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

ObjectProvider<ArcJetLTEPhysicalModel<MultiScalarTerm<EulerTerm>, NSTerm, DIM_2D>, 
	       PhysicalModelImpl, ArcJetModule, 1>
arcJetLTE2DProvider("ArcJetLTE2D");

ObjectProvider<ArcJetLTEPhysicalModel<MultiScalarTerm<EulerTerm>, NSTerm, DIM_3D>, 
	       PhysicalModelImpl, ArcJetModule, 1>
arcJetLTE3DProvider("ArcJetLTE3D");

ObjectProvider<ArcJetLTEPhysicalModel<MultiScalarTerm<EulerTerm>, NSTerm, DIM_2D>, 
	       PhysicalModelImpl, ArcJetModule, 1>
arcJetPG2DProvider("ArcJetPG2D");

ObjectProvider<ArcJetLTEPhysicalModel<MultiScalarTerm<EulerTerm>, NSTerm, DIM_3D>, 
	       PhysicalModelImpl, ArcJetModule, 1>
arcJetPG3DProvider("ArcJetPG3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
