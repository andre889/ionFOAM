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

#include "CompactListList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
template<class Container>
Foam::CompactListList<T>::CompactListList(const List<Container>& ll)
:
    UCompactListList<T>(),
    offsets_(ll.size() + 1),
    m_()
{
    label sumSize = 0;
    offsets_[0] = 0;
    forAll(ll, i)
    {
        sumSize += ll[i].size();
        offsets_[i+1] = sumSize;
    }

    m_.setSize(sumSize);

    label k = 0;
    forAll(ll, i)
    {
        forAll(ll[i], j)
        {
            m_[k++] = ll[i][j];
        }
    }

    UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    const labelUList& rowSizes,
    const T& t
)
:
    UCompactListList<T>(),
    offsets_(rowSizes.size() + 1),
    m_()
{
    label sumSize = 0;
    offsets_[0] = 0;
    forAll(rowSizes, i)
    {
        sumSize += rowSizes[i];
        offsets_[i+1] = sumSize;
    }

    m_.setSize(sumSize, t);

    UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::CompactListList<T>::setSize(const label mRows)
{
    if (mRows == 0)
    {
        clear();
    }
    if (mRows < this->size())
    {
        offsets_.setSize(mRows + 1);
        m_.setSize(offsets_[mRows]);

        UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));
    }
    else if (mRows > this->size())
    {
        FatalErrorInFunction
            << "Cannot be used to extend the list from " << offsets_.size()
            << " to " << mRows << nl
            << "    Please use one of the other setSize member functions"
            << abort(FatalError);
    }
}


template<class T>
void Foam::CompactListList<T>::setSize
(
    const label mRows,
    const label nData
)
{
    offsets_.setSize(mRows + 1);
    m_.setSize(nData);

    UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));
}


template<class T>
void Foam::CompactListList<T>::setSize
(
    const label mRows,
    const label nData,
    const T& t
)
{
    offsets_.setSize(mRows + 1);
    m_.setSize(nData, t);

    UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));
}


template<class T>
void Foam::CompactListList<T>::setSize(const labelUList& rowSizes)
{
    offsets_.setSize(rowSizes.size() + 1);

    label sumSize = 0;
    offsets_[0] = 0;
    forAll(rowSizes, i)
    {
        sumSize += rowSizes[i];
        offsets_[i+1] = sumSize;
    }

    m_.setSize(sumSize);

    UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));
}


template<class T>
void Foam::CompactListList<T>::clear()
{
    offsets_ = List<label>(1, 0);
    m_.clear();

    UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));
}


template<class T>
void Foam::CompactListList<T>::transfer(CompactListList<T>& a)
{
    offsets_.transfer(a.offsets_);
    m_.transfer(a.m_);

    UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));

    a.clear();
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CompactListListIO.C"

// ************************************************************************* //
