Info << "Reading velocity field\n"
     << endl;

volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        meshFluid,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    meshFluid);


surfaceScalarField phi(
    IOobject(
        "phi",
        runTime.timeName(),
        meshFluid,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    linearInterpolate(U) & meshFluid.Sf());

volScalarField p(
    IOobject(
        "p",
        runTime.timeName(),
        meshFluid,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    meshFluid);

volVectorField Fe(
IOobject(
        "Fe",
        runTime.timeName(),
        meshFluid,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
        meshFluid,
        dimensionedVector("0", p.dimensions()/dimLength,vector(0,0,0)),
        zeroGradientFvPatchScalarField::typeName
        );

//Hmmm BCS?
surfaceScalarField* phiGlobalPtr;
volVectorField UGlobal(
    IOobject(
        "Uglobal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedVector("0", U.dimensions(), vector(0, 0, 0)),
    zeroGradientFvPatchScalarField::typeName);
{

    /*PtrList<fvPatchVectorField> bcList;
    label extrBoundaryIndex = mesh.boundaryMesh().size();
    forAll(meshFluid.boundaryMesh(),fluidBoundaryi)
    {
        if(boundaryRegionAddressList[fluidBoundaryi] > -1)
        {
            bcList.set(fluidBoundaryi, new fvPatchVectorField(U.boundaryFieldRef()[fluidBoundaryi]));
        }
        else
        {

        bcList.set(extrBoundaryIndex, new fvPatchVectorField(zeroGradientFvPatchField));
        extrBoundaryIndex++;
        }
    }*/

    
/*
forAll(meshFluid.boundaryMesh(),fluidBoundaryi)
            {
                fvPatchVectorField& bci = U.boundaryFieldRef()[fluidBoundaryi]; 
                label startIndex = bci.patch().start();
                label sizeOfPatch = bci.patch().size(); 
                //If not -1, then boundary exists in global mesh. 
                //Otherwise, fluid boundary data needs to be mapped to internal faces of phiGlobal.
                //Info << "start index " << startIndex << " size " << sizeOfPatch << " for bc " << fluidBoundaryi << endl;

                if(boundaryRegionAddressList[fluidBoundaryi] > -1)
                {
                       fvPatchVectorField& bcGlobali = UGlobal.boundaryFieldRef()[boundaryRegionAddressList[fluidBoundaryi]];
                       bcGlobali = bci;
                        Info << "Copying Boundary: " << fluidBoundaryi << " of Type: "<< bci.typeName;
                         Info << "To Boundary: " << boundaryRegionAddressList[fluidBoundaryi] << " of Type: "<< bcGlobali.typeName;

                }
            }*/
        



    phiGlobalPtr = new surfaceScalarField(
    IOobject(
        "phiGlobal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
        linearInterpolate(UGlobal) & mesh.Sf());
}

surfaceScalarField& phiGlobal = *phiGlobalPtr;