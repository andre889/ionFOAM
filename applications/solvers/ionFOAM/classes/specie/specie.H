//  Date:  05/05/2023
//  Revision: 1.0
//  Author: Austin J. Andrews
//  Email: andre889@umn.edu
//  Affiliation: University of Minnesota Twin Cities, Department of Mechanical Engineering, Hogan Lab

//  Description: This class allows for tracking of specie variables and their associated properties such as charge, diffusion coefficient and the single ion partition coefficient.
// **********************************************

#ifndef specie_H
#define specie_H

#include "fvCFD.H"


namespace Foam
{
    class specie
    {
    private:
        volScalarField s_; //Field data for specie i
        const dictionary& dict_;

    public:
        //Ion properties
        scalar z; // charge state
        scalar Daq; // diffusion coefficient in water (reference solution)
        scalar Dmem; // diffusion coefficient in membrane (reference solution)
        scalar K; // single ion partition coefficient 
        scalar C0membrane; // Initial concentration in membrane
        scalar C0aq; // initial concentration in water (reference solution)

        specie(IOobject IOobj, const fvMesh& mesh,const dictionary &subDict,bool nothing);
        specie(IOobject IOobj, const fvMesh &mesh, const dictionary &subDict);

        volScalarField &s() { return s_; };

        ~specie();
    };
    
}

#endif