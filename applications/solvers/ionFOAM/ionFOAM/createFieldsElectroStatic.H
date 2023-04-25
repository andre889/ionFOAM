
Info << "Create Electrostatics\n";
//**********************************
//Electrostatic Properties
//*********************************
volScalarField epsilon(
     IOobject(
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    constant::electromagnetic::epsilon0,
    zeroGradientFvPatchScalarField::typeName);

//e = e*erAq(1-alpha) + (e*erMem)(alpha)
epsilon = epsilon*liquid.epsilon.value()*(1.0-alpha) + epsilon*membrane.epsilon.value()*(alpha);


//electroStatic poissonFields(z,epsilon,mesh);
volScalarField q(
    IOobject(
        "q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimMoles/dimVolume*constant::physicoChemical::F.dimensions(), 0),
    zeroGradientFvPatchScalarField::typeName);

    volScalarField voltage(
    IOobject(
        "voltage",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

label voltageRefCell = 0;
scalar voltageRefValue = 0;
setRefCell(voltage, piso.dict(), voltageRefCell, voltageRefValue);

volVectorField E(IOobject(
                     "E",
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::NO_WRITE),
                 -fvc::grad(voltage));


Info << "Done with Create Fields\n";
