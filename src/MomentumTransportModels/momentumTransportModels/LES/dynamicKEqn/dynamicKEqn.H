/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
     Foam::LESModels::dynamicKEqn

Description
    Dynamic one equation eddy-viscosity model

    Eddy viscosity SGS model using a modeled balance equation to simulate
    the behaviour of k in which a dynamic procedure is applied to evaluate the
    coefficients.

    Reference:
    \verbatim
        Kim, W and Menon, S. (1995).
        A new dynamic one-equation subgrid-scale model for
        large eddy simulation.
        In 33rd Aerospace Sciences Meeting and Exhibit, Reno, NV, 1995.
    \endverbatim

    There are no default model coefficients but the filter used for KK must be
    supplied, e.g.
    \verbatim
        dynamicKEqnCoeffs
        {
            filter simple;
        }
    \endverbatim

SourceFiles
    dynamicKEqn.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicKEqn_H
#define dynamicKEqn_H

#include "LESeddyViscosity.H"
#include "simpleFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class dynamicKEqn Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class dynamicKEqn
:
    public LESeddyViscosity<BasicMomentumTransportModel>
{
protected:

    // Protected data

        // Fields

            volScalarField k_;


        // Filters

            simpleFilter simpleFilter_;
            autoPtr<LESfilter> filterPtr_;
            LESfilter& filter_;


    // Protected Member Functions

        //- Return the KK field obtained by filtering the velocity field U
        tmp<volScalarField> KK() const;

        //- Calculate Ck by filtering the velocity field U
        volScalarField Ck
        (
            const volSymmTensorField& D,
            const volScalarField& KK
        ) const;

        //- Calculate Ce by filtering the velocity field U
        volScalarField Ce
        (
            const volSymmTensorField& D,
            const volScalarField& KK
        ) const;

        volScalarField Ce() const;

        //- Update sub-grid eddy-viscosity
        void correctNut
        (
            const volSymmTensorField& D,
            const volScalarField& KK
        );

        virtual void correctNut();

        virtual tmp<fvScalarMatrix> kSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("dynamicKEqn");


    // Constructors

        //- Construct from components
        dynamicKEqn
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        dynamicKEqn(const dynamicKEqn&) = delete;


    //- Destructor
    virtual ~dynamicKEqn()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                this->nut_ + this->nu()
            );
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicKEqn&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicKEqn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
