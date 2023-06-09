/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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
    Foam::LESModels::continuousGasKEqn

Description
    One-equation SGS model for the gas-phase in a two-phase system
    supporting phase-inversion.

    In the limit that the gas-phase fraction approaches zero a contribution from
    the other phase is blended into the k-equation up to the phase-fraction of
    alphaInversion at which point phase-inversion is considered to have occurred
    and the model reverts to the pure single-phase form.

    This model is unpublished and is provided as a stable numerical framework
    on which a more physical model may be built.

    The default model coefficients are
    \verbatim
        continuousKEqnCoeffs
        {
            Ck              0.094;
            Ce              1.048;
            alphaInversion  0.7;
        }
    \endverbatim

SourceFiles
    continuousGasKEqn.C

\*---------------------------------------------------------------------------*/

#ifndef continuousGasKEqn_H
#define continuousGasKEqn_H

#include "kEqn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                      Class continuousGasKEqn Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class continuousGasKEqn
:
    public kEqn<BasicMomentumTransportModel>
{
    // Private Data

        mutable const momentumTransportModel *liquidTurbulencePtr_;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar alphaInversion_;


    // Protected Member Functions

        //- Return the turbulence model for the liquid phase
        const momentumTransportModel& liquidTurbulence() const;

        tmp<volScalarField> phaseTransferCoeff() const;
        virtual tmp<fvScalarMatrix> kSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("continuousGasKEqn");


    // Constructors

        //- Construct from components
        continuousGasKEqn
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
        continuousGasKEqn(const continuousGasKEqn&) = delete;


    //- Destructor
    virtual ~continuousGasKEqn()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const continuousGasKEqn&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "continuousGasKEqn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
