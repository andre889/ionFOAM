//***********************************************
//  Date:  05/05/2023
//  Revision: 1.0
//  Author: Austin J. Andrews
//  Email: andre889@umn.edu
//  Affiliation: University of Minnesota Twin Cities, Department of Mechanical Engineering, Hogan Lab

//  Description: See myddtScheme.H
// **********************************************

#include "myddtScheme.H"

namespace Foam
{
    //Set upsSelection table, but do not include myddtScheme as a selection.
      defineTypeNameAndDebug(myddtScheme, 0);
      defineRunTimeSelectionTable(myddtScheme, Normal);      

//Constructor for all myddtSchemes
myddtScheme::myddtScheme
    (
        const fvMesh& mesh,
        const dictionary& dict
    )
    :
    mesh_(mesh),
    dict_(dict)
    {}
       
//Create a new ddt scheme from a list of implemented schemes.
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

