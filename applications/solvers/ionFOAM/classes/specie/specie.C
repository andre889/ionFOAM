#include "specie.H"

namespace Foam
{
    
        specie::specie(IOobject IOobj, const fvMesh& mesh,const dictionary &subDict,bool nothing):
        s_(IOobj,mesh,subDict),
        dict_(subDict)
        {
        Info << "Gathering Properties for ion: " << s_.name() << endl;

        z = readScalar(dict_.lookup("z"));
        Dmem = readScalar(dict_.lookup("Dmem"));
        Daq = readScalar(dict_.lookup("Daq"));
        K = readScalar(dict_.lookup("K"));
        C0membrane = readScalar(dict_.lookup("C0membrane"));
        C0aq = readScalar(dict_.lookup("C0aq"));

        Info << "K" << K << endl;
        
        //Info << "Reading BoundaryField for " << s_.name() << endl;
        //s_.boundaryFieldRef().readField(s_, dict_.subDict("boundaryField"));

        };

        specie::specie(IOobject IOobj, const fvMesh &mesh, const dictionary &subDict) : s_(IOobj, mesh),
                                                                                                      dict_(subDict)
        {
                Info << "Gathering Properties for ion: " << s_.name() << endl;

                z = readScalar(dict_.lookup("z"));
                Dmem = readScalar(dict_.lookup("Dmem"));
                Daq = readScalar(dict_.lookup("Daq"));
                K = readScalar(dict_.lookup("K"));
                C0membrane = readScalar(dict_.lookup("C0membrane"));
                C0aq = readScalar(dict_.lookup("C0aq"));

                Info << "K" << K << endl;

                //Info << "Reading BoundaryField for " << s_.name() << endl;
                //s_.boundaryFieldRef().readField(s_, dict_.subDict("boundaryField"));
        };

        specie::~specie(){};
    
    
}
