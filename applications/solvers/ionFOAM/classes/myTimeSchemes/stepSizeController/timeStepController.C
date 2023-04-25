#include "timeStepController.H"
//fk


namespace Foam
{
    //Set upsSelection table, but do not include
    //aerosolDistribution as a selection.
      defineTypeNameAndDebug(timeStepController, 0);
      defineRunTimeSelectionTable(timeStepController, Normal);      

timeStepController::timeStepController
    (
        const fvMesh& mesh,
        const dictionary& dict,
        scalar p
    )
    :
    mesh_(mesh),
    dict_(dict),
    runTime_(const_cast<Time &>(mesh.time()))
    {
        phat = p-1;
        epsilon = readScalar(dict.lookup("error"));
        minDt = readScalar(dict.lookup("minDt"));
        maxRelativeError = dict.lookupOrDefault<scalar>("maxRelativeError", 10);
        dtLimitFactor = dict.lookupOrDefault<scalar>("dtLimitFactor", 5);
        maxDt = dict.lookupOrDefault<scalar>("maxDt", 1e3);
        K = dict.lookupOrDefault<scalar>("safteyFactor", 0.95);
        K = min(K, 1.0);
        K = max(K, 0.5);
    }
       
             
autoPtr<timeStepController> timeStepController::New(
        const fvMesh& mesh,
        const dictionary& dict,
        scalar p
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        dict.lookup("timeStepController")
    );

    Info<< "Selecting stepController " << modelType << endl;

    typename NormalConstructorTable::iterator cstrIter
        = NormalConstructorTablePtr_->find(modelType);
   
    if (cstrIter == NormalConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown ddtScheme type "
            << modelType << nl << nl
            << "Valid ddtScheme types:" << endl
            << NormalConstructorTablePtr_->sortedToc()
            << "\n░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ \n░░░░░░░░░░░░░░░░▄▄███▄▄▄░▄▄██▄░░░░░░░ \n░░░░░░░░░██▀███████████████▀▀▄█░░░░░░ \n░░░░░░░░█▀▄▀▀▄██████████████▄░░█░░░░░ \n░░░░░░░█▀▀░▄██████████████▄█▀░░▀▄░░░░ \n░░░░░▄▀░░░▀▀▄████████████████▄░░░█░░░ \n░░░░░▀░░░░▄███▀░░███▄████░████░░░░▀▄░ \n░░░▄▀░░░░▄████░░▀▀░▀░░░░░░██░▀▄░░░░▀▄ \n░▄▀░░░░░▄▀▀██▀░░░░░▄░░▀▄░░██░░░▀▄░░░░ \n█░░░░░█▀░░░██▄░░░░░▀▀█▀░░░█░░░░░░█░░░ \n█░░░▄▀░░░░░░██░░░░░▀██▀░░█▀▄░░░░░░▀▀▀ \n▀▀▀▀░▄▄▄▄▄▄▀▀░█░░░░░░░░░▄█░░█▀▀▀▀▀█░░ \n░░░░█░░░▀▀░░░░░░▀▄░░░▄▄██░░░█░░░░░▀▄░ \n░░░░█░░░░░░░░░░░░█▄▀▀▀▀▀█░░░█░░░░░░█░ \n░░░░▀░░░░░░░░░░░░░▀░░░░▀░░░░▀░░░░░░░░ \n░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ \n"
            << exit(FatalError);
    }
    
    return autoPtr<timeStepController> (cstrIter()(mesh,dict,p));
}

timeStepController::~timeStepController(){}


}

