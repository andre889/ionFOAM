/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "GidaspowSchillerNaumann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(GidaspowSchillerNaumann, 0);
    addToRunTimeSelectionTable(dragModel, GidaspowSchillerNaumann, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::GidaspowSchillerNaumann::GidaspowSchillerNaumann
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject),
    residualRe_("residualRe", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::GidaspowSchillerNaumann::~GidaspowSchillerNaumann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::GidaspowSchillerNaumann::CdRe() const
{
    const volScalarField alpha2
    (
        max(1 - interface_.dispersed(), interface_.continuous().residualAlpha())
    );

    const volScalarField Re(alpha2*interface_.Re());

    const volScalarField CdsRe
    (
        neg(Re - 1000)*24*(1.0 + 0.15*pow(Re, 0.687))/alpha2
      + pos0(Re - 1000)*0.44*max(Re, residualRe_)
    );

    return
        CdsRe
       *pow(alpha2, -2.65)
       *max(interface_.continuous(), interface_.continuous().residualAlpha());
}


// ************************************************************************* //
