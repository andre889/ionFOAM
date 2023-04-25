#include "PPIDStepController.H"
//fk


namespace Foam
{

    defineTypeNameAndDebug(PPIDStepController, 0);
    addToRunTimeSelectionTable(timeStepController, PPIDStepController, Normal);

PPIDStepController::PPIDStepController
    (
        const fvMesh& mesh,
        const dictionary& dict,
        scalar p
    )
    :
    timeStepController(mesh,dict,p)
    {
        Info << "hello from PC constructor \n";

    }
       


PPIDStepController::~PPIDStepController(){}


}

