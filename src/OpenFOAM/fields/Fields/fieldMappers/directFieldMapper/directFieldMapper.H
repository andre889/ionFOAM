/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2022 OpenFOAM Foundation
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
    Foam::directFieldMapper

Description
    Direct field mapper

\*---------------------------------------------------------------------------*/

#ifndef directFieldMapper_H
#define directFieldMapper_H

#include "fieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class directFieldMapper
:
    virtual public fieldMapper
{
    // Private Member Functions

        template<class Type>
        void map(Field<Type>& f, const Field<Type>& mapF) const;

        template<class Type>
        tmp<Field<Type>> map(const Field<Type>& f) const;


public:

    // Constructors

        //- Null constructor
        directFieldMapper()
        {}


    //- Destructor
    virtual ~directFieldMapper()
    {}


    // Member Functions

        //- Are there unmapped values? I.e. do all size() elements get
        //  get value
        virtual bool hasUnmapped() const = 0;

        //- Access to the direct map addressing
        virtual const labelUList& addressing() const = 0;


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "directFieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
