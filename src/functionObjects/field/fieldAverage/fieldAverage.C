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

\*---------------------------------------------------------------------------*/

#include "fieldAverage.H"
#include "fieldAverageItem.H"
#include "timeIOdictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldAverage, 0);
    addToRunTimeSelectionTable(functionObject, fieldAverage, dictionary);
}
}


template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldAverage::baseType,
    2
>::names[] = { "iteration", "time"};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldAverage::baseType,
    2
> Foam::functionObjects::fieldAverage::baseTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldAverage::resetFields()
{
    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            if (obr_.found(faItems_[i].meanFieldName()))
            {
                obr_.checkOut(*obr_[faItems_[i].meanFieldName()]);
            }
        }

        if (faItems_[i].prime2Mean())
        {
            if (obr_.found(faItems_[i].prime2MeanFieldName()))
            {
                obr_.checkOut(*obr_[faItems_[i].prime2MeanFieldName()]);
            }
        }
    }
}


void Foam::functionObjects::fieldAverage::initialise()
{
    if (!totalIter_.size())
    {
        totalIter_.setSize(faItems_.size(), 1);
    }

    if (!totalTime_.size())
    {
        totalTime_.setSize(faItems_.size(), obr_.time().deltaTValue());
    }
    else
    {
        // Check if totalTime_ has been set otherwise initialise
        forAll(totalTime_, fieldi)
        {
            if (totalTime_[fieldi] < 0)
            {
                totalTime_[fieldi] = obr_.time().deltaTValue();
            }
        }
    }

    resetFields();

    // Add mean fields to the field lists
    forAll(faItems_, fieldi)
    {
        addMeanField<scalar>(fieldi);
        addMeanField<vector>(fieldi);
        addMeanField<sphericalTensor>(fieldi);
        addMeanField<symmTensor>(fieldi);
        addMeanField<tensor>(fieldi);
    }

    // Add prime-squared mean fields to the field lists
    forAll(faItems_, fieldi)
    {
        addPrime2MeanField<scalar, scalar>(fieldi);
        addPrime2MeanField<vector, symmTensor>(fieldi);
    }

    // ensure first averaging works unconditionally
    prevTimeIndex_ = -1;

    initialised_ = true;
}


void Foam::functionObjects::fieldAverage::restart()
{
    Log << "    Restarting averaging at time " << obr_.time().timeName()
        << nl << endl;

    totalIter_.clear();
    totalTime_.clear();

    initialise();
}


void Foam::functionObjects::fieldAverage::calcAverages()
{
    if (!initialised_)
    {
        initialise();
    }

    const label currentTimeIndex = obr_.time().timeIndex();
    const scalar currentTime = obr_.time().value();

    if (prevTimeIndex_ == currentTimeIndex)
    {
        return;
    }
    else
    {
        prevTimeIndex_ = currentTimeIndex;
    }

    if (periodicRestart_ && currentTime > restartPeriod_*periodIndex_)
    {
        restart();
        periodIndex_++;
    }

    Log << type() << " " << name() << nl
        << "    Calculating averages" << nl;

    addMeanSqrToPrime2Mean<scalar, scalar>();
    addMeanSqrToPrime2Mean<vector, symmTensor>();

    calculateMeanFields<scalar>();
    calculateMeanFields<vector>();
    calculateMeanFields<sphericalTensor>();
    calculateMeanFields<symmTensor>();
    calculateMeanFields<tensor>();

    calculatePrime2MeanFields<scalar, scalar>();
    calculatePrime2MeanFields<vector, symmTensor>();

    forAll(faItems_, fieldi)
    {
        totalIter_[fieldi]++;
        totalTime_[fieldi] += obr_.time().deltaTValue();
    }

    Log << endl;
}


void Foam::functionObjects::fieldAverage::writeAverages() const
{
    Log << "    Writing average fields" << endl;

    writeFields<scalar>();
    writeFields<vector>();
    writeFields<sphericalTensor>();
    writeFields<symmTensor>();
    writeFields<tensor>();

    Log << endl;
}


void Foam::functionObjects::fieldAverage::writeAveragingProperties() const
{
    timeIOdictionary propsDict
    (
        IOobject
        (
            name() + "Properties",
            obr_.time().timeName(),
            "uniform",
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(faItems_, fieldi)
    {
        const word& fieldName = faItems_[fieldi].fieldName();
        propsDict.add(fieldName, dictionary());
        propsDict.subDict(fieldName).add("totalIter", totalIter_[fieldi]);
        propsDict.subDict(fieldName).add("totalTime", totalTime_[fieldi]);
    }

    propsDict.regIOobject::write();

    Log << endl;
}


void Foam::functionObjects::fieldAverage::readAveragingProperties()
{
    if ((restartOnRestart_ || restartOnOutput_) && log)
    {
        Info<< "    Starting averaging at time " << obr_.time().timeName()
            << nl;
    }
    else
    {
        typeIOobject<timeIOdictionary> propsDictHeader
        (
            name() + "Properties",
            obr_.time().timeName(obr_.time().startTime().value()),
            "uniform",
            obr_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        );

        if (!propsDictHeader.headerOk())
        {
            Log << "    Starting averaging at time "
                << obr_.time().timeName() << nl;

            return;
        }

        timeIOdictionary propsDict(propsDictHeader);

        Log << "    Restarting averaging for fields:" << nl;

        totalIter_.setSize(faItems_.size(), 1);

        // Initialise totalTime with negative values
        // to indicate that it has not been set
        totalTime_.setSize(faItems_.size(), -1);

        forAll(faItems_, fieldi)
        {
            const word& fieldName = faItems_[fieldi].fieldName();
            if (propsDict.found(fieldName))
            {
                dictionary fieldDict(propsDict.subDict(fieldName));

                totalIter_[fieldi] = fieldDict.lookup<label>("totalIter");
                totalTime_[fieldi] = fieldDict.lookup<scalar>("totalTime");

                Log << "        " << fieldName
                    << " iters = " << totalIter_[fieldi]
                    << " time = " << totalTime_[fieldi] << nl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverage::fieldAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    prevTimeIndex_(-1),
    restartOnRestart_(false),
    restartOnOutput_(false),
    periodicRestart_(false),
    restartPeriod_(great),
    initialised_(false),
    base_(baseType::iter),
    window_(-1.0),
    windowName_(""),
    faItems_(),
    totalIter_(),
    totalTime_(),
    periodIndex_(1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverage::~fieldAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldAverage::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    initialised_ = false;

    Log << type() << " " << name() << ":" << nl;

    dict.readIfPresent("restartOnRestart", restartOnRestart_);
    dict.readIfPresent("restartOnOutput", restartOnOutput_);
    dict.readIfPresent("periodicRestart", periodicRestart_);

    mean_ = dict.lookupOrDefault<Switch>("mean", true);
    prime2Mean_ = dict.lookupOrDefault<Switch>("prime2Mean", false);
    base_ = baseTypeNames_
    [
        dict.lookupOrDefault<word>("base", "time")
    ];
    window_ = dict.lookupOrDefault<scalar>("window", -1);
    windowName_ = dict.lookupOrDefault<word>("windowName", "");

    faItems_ = PtrList<fieldAverageItem>
    (
        dict.lookup("fields"),
        fieldAverageItem::iNew(*this)
    );

    if (periodicRestart_)
    {
        dict.lookup("restartPeriod") >> restartPeriod_;
    }

    readAveragingProperties();

    Log << endl;

    return true;
}


Foam::wordList Foam::functionObjects::fieldAverage::fields() const
{
    wordList fields(faItems_.size());

    forAll(faItems_, fieldi)
    {
        fields[fieldi] = faItems_[fieldi].fieldName();
    }

    return fields;
}


bool Foam::functionObjects::fieldAverage::execute()
{
    calcAverages();

    return true;
}


bool Foam::functionObjects::fieldAverage::write()
{
    writeAverages();
    writeAveragingProperties();

    if (restartOnOutput_)
    {
        restart();
    }

    return true;
}


// ************************************************************************* //
