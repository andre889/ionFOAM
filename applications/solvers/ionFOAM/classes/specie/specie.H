#ifndef specie_H
#define specie_H

#include "fvCFD.H"


namespace Foam
{
    class specie
    {
    private:
        volScalarField s_; //main specie class
        PtrList<volScalarField> Ki_; //runge kutta List
        volScalarField* oldTimePtr; 
        const dictionary& dict_;

        //Later add diffusion coeff

    public:
    //Ion properties
        scalar z; // charge state
         scalar Daq; // diffusion coefficient in water (reference solution)
        scalar Dmem; // diffusion coefficient in membrane (reference solution)
        scalar K; // partition coefficient 
        scalar C0membrane; // Initial concentration in membrane
        scalar C0aq; // intiial concentration in water (reference solution)

        specie(IOobject IOobj, const fvMesh& mesh,const dictionary &subDict,bool nothing);
        specie(IOobject IOobj, const fvMesh &mesh, const dictionary &subDict);

        volScalarField &s() { return s_; };

        ~specie();
    };
    
}

#endif