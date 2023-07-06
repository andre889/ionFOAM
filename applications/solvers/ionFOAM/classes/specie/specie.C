//  Date:  05/05/2023
//  Revision: 1.0
//  Author: Austin J. Andrews
//  Email: andre889@umn.edu
//  Affiliation: University of Minnesota Twin Cities, Department of Mechanical Engineering, Hogan Lab

//  Description: See specie.H
// **********************************************

#include "specie.H"

namespace Foam
{

        specie::specie(IOobject IOobj, const fvMesh &mesh, const dictionary &subDict, bool nothing) : s_(IOobj, mesh, subDict),
                                                                                                      dict_(subDict)
        {
                Info << "Gathering Properties for ion: " << s_.name() << endl;

                z = readScalar(dict_.lookup("z"));
                Dmem = readScalar(dict_.lookup("Dmem"));
                Daq = readScalar(dict_.lookup("Daq"));
                K = readScalar(dict_.lookup("K"));
                C0membrane = readScalar(dict_.lookup("C0membrane"));
                C0aq = readScalar(dict_.lookup("C0aq"));
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
        };

        specie::~specie(){};

}
