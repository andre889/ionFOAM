#include "myddtScheme.H"
//fk


namespace Foam
{
    //Set upsSelection table, but do not include
    //aerosolDistribution as a selection.
      defineTypeNameAndDebug(myddtScheme, 0);
      defineRunTimeSelectionTable(myddtScheme, Normal);      

myddtScheme::myddtScheme
    (
        const fvMesh& mesh,
        const dictionary& dict
    )
    :
    mesh_(mesh),
    dict_(dict)
    {
        //Constructor for all aerosolDistributions
        
    }
       
             
autoPtr<myddtScheme> myddtScheme::New(
        const fvMesh& mesh,
        const dictionary& dict
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        dict.lookup("myddtScheme")
    );

    Info<< "Selecting ddtScheme " << modelType << endl;

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
    
    return autoPtr<myddtScheme> (cstrIter()(mesh,dict));
}

myddtScheme::~myddtScheme(){}


}

