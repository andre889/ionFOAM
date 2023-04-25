
word fluidRegionName = "fluid";
word membraneRegionName = "membrane";
word totalRegionName =  fvMesh::defaultRegion;//totalRegionName



Info << "Create mesh for time = "
    << runTime.timeName() << nl << endl;

bool multiRegion = true;
regionProperties rp(runTime);

Info << rp << endl;

Info << totalRegionName << endl;


fvMesh* meshFluidPtr;
fvMesh* meshPtr; 

//default mesh
fvMesh meshReg0(IOobject(
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ));

meshPtr = &meshReg0;

if(rp.size() == 1)
{
    Info << " ONLY 1 region! " << endl;
    fluidRegionName = fvMesh::defaultRegion;
    multiRegion = false;
    meshFluidPtr = &meshReg0;
}
else if(rp.size() == 2)
{
    
meshFluidPtr = new fvMesh(
    IOobject(
        fluidRegionName,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ));
}

else if(rp.size() > 2)
{
    FatalError  << " More than 2 regions not yet supported " << exit(FatalError);
}


fvMesh& meshFluid = *meshFluidPtr; 
fvMesh& mesh = *meshPtr;

Info << "Done creating mesh!\n";
Info << "Default Mesh has: " << mesh.C().size() << " cells\n";

//determin if multiple zones are pressent. If not then 



Info << "Fluid Mesh has: " << meshFluid.C().size() << " cells\n";

Info << "Creating addressing for fluid to global mesh regions...\n";


  std::vector<label> cellAddressListFluid;
  std::vector<label> boundaryRegionAddressList;
  std::vector<label> faceAddressList;
      
  
      
if(multiRegion)
{
    
fvMesh& meshRef = meshFluid;
fileName directoryOfMesh = meshRef.meshDir();
meshRef.owner();
Info << "Directory of fluid Mesh: " << directoryOfMesh << endl;
fileName newDir = runTime.time().constant() +"/" + directoryOfMesh;
Info << "Directory of faceAddress: " << newDir << endl;

 labelIOList faceAddressing1(
        IOobject(
            "polyMesh/faceRegionAddressing",
            runTime.time().constant(),
            meshFluid,
            IOobject::MUST_READ,
            IOobject::NO_WRITE));

    labelIOList boundaryRegionAddressing(
        IOobject(
            "polyMesh/boundaryRegionAddressing",
            runTime.time().constant(),
            meshFluid,
            IOobject::MUST_READ,
            IOobject::NO_WRITE));

    labelIOList cellAddressing(
      IOobject(
                "polyMesh/cellRegionAddressing",
                runTime.time().constant(),
                meshFluid,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
     );
     
     cellAddressListFluid.resize(cellAddressing.size());
     boundaryRegionAddressList.resize(boundaryRegionAddressing.size());
    faceAddressList.resize(faceAddressing1.size());

    for (label i = 0; i < faceAddressing1.size(); ++i)
    {
        faceAddressList[i] = faceAddressing1[i];
    }
     for (label i = 0; i < boundaryRegionAddressing.size(); ++i)
    {
        boundaryRegionAddressList[i] = boundaryRegionAddressing[i];
    }
     for (label i = 0; i < cellAddressing.size(); ++i)
    {
        cellAddressListFluid[i] = cellAddressing[i];
    }

    //for some reason parallel cases get the boundaryRegionAddressList incorrect. To fix this I will manually reconstruct it. 
        boundaryRegionAddressList.clear();
        boundaryRegionAddressList.resize(meshFluid.boundaryMesh().names().size());
        
        forAll(meshFluid.boundaryMesh().names(),fluidBoundaryi)
        {
            boundaryRegionAddressList[fluidBoundaryi] = -1;
          word boundaryNamei = meshFluid.boundaryMesh().names()[fluidBoundaryi];
          forAll(mesh.boundaryMesh().names(),fluidBoundaryj)
          {
            word boundaryNamej = mesh.boundaryMesh().names()[fluidBoundaryj];
             if( boundaryNamej== boundaryNamei)
             {
                boundaryRegionAddressList[fluidBoundaryi] = fluidBoundaryj;
                break;
             }
          }
        }
}
else //test for parallel?
{ 
    if(Pstream::parRun())
    {
        Info << "Warning Not tested for Parallel! LINE 152\n ";
    }

    cellAddressListFluid.resize(mesh.C().size());
    for (label i = 0; i < mesh.C().size(); ++i)
    {
        cellAddressListFluid[i] = i;
    }
}
  
   
//Create Mesh with two regions 
//One for the whole mesh to solve ion transport and poissons equation
//One for just the fluid region to solve for the NS equations with coupled body force term!
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
        IOobject::NO_WRITE),
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
        IOobject::NO_WRITE),
        meshFluid,
        dimensionedVector("0", p.dimensions()/dimLength,vector(0,0,0)),
        zeroGradientFvPatchScalarField::typeName
        );

//Hmmm BCS?
surfaceScalarField* phiGlobalPtr;
//{ delete for debug

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

    volVectorField UGlobal (
                    IOobject(
                            "Uglobal",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE),
                            mesh,
                            dimensionedVector("0", U.dimensions(),vector(0,0,0)),
                            zeroGradientFvPatchScalarField::typeName);
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
        IOobject::NO_WRITE),
        linearInterpolate(UGlobal) & mesh.Sf());
//} delete for debug

surfaceScalarField& phiGlobal = *phiGlobalPtr;

//Create Solution Control
pisoControl piso(mesh);
//For now just used for dictionary
const dictionary& pisoDictRef = piso.dict();


label pRefCell = 0;
scalar pRefValue = 0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);


meshFluid.setFluxRequired(p.name());


autoPtr<viscosityModel> viscosity(viscosityModel::New(meshFluid));

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, viscosity)
);

//*****************************
//Gathering Mat'l properties
//*****************************
Info << "Reading materialProperties\n"
     << endl;

//Dictionary Definition
IOdictionary materialProperties(
    IOobject(
        "materialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE));



const dictionary &membraneSubDict = materialProperties.subDict("Membrane");
material membrane(membraneSubDict);


const dictionary &liqiudSubDict = materialProperties.subDict("Liquid");
material liquid(liqiudSubDict);


Info << "Reading transportProperties\n"
     << endl;

//Dictionary Definition
IOdictionary transportProperties(
    IOobject(
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE));

dimensionedScalar rhof(transportProperties.lookup("rhof"));
dimensionedScalar T(transportProperties.lookup("T"));

bool electroPhoretic = readBool(transportProperties.lookup("electroPhoretic"));
bool diffusion  = readBool(transportProperties.lookup("diffusion"));
bool chemcialPotential = readBool(transportProperties.lookup("chemcialPotential"));


Info << "Reading ionProperties\n"
     << endl;

//Dictionary Definition
IOdictionary ionProperties(
    IOobject(
        "ionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE));

//List of ions to use
wordList ionNames(
    ionProperties.lookup("ions"));
label numIons = ionNames.size();

//Generate Lists for each species
PtrList<specie> Ci(numIons); //Concentration

//https://wiki.openfoam.com/Programming2 set IC easily!
// readFields(dict)

forAll(Ci, n)
{   
    word currentIonName = ionNames[n];

    //If Ion exists in main director read from this file. Otherwise Create from transport property
    //Force overwrite?
    IOobjectList objects(mesh, runTime.timeName());
    wordList objectNames = objects.names();
    Info << "Trying to find Ion named: " << currentIonName << endl;
    //wordList::iterator namePtr = objectNames.find(currentIonName);
    label index = -1;
    forAll(objectNames,oi)
    {
        if (objectNames[oi] == currentIonName) //compare
        {
            index = oi;
            break;
        }
    }

    bool ionFound = (index != -1);
    // bool ionFound = objectNames
    Info << "Ion Found ? " << ionFound << endl;

    const dictionary &subDict = ionProperties.subDict(currentIonName);

    if (ionFound)
    {
        Ci.set(n, new specie(
                      IOobject(
                          ionNames[n],
                          runTime.timeName(),
                          mesh,
                          IOobject::MUST_READ,
                          IOobject::AUTO_WRITE),
                            mesh,subDict));
    }
    else
    {
        Ci.set(n, new specie(
                    IOobject(
                        ionNames[n],
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE),
                        mesh,subDict
                    ,true));

        //possible error with dimset
        Ci[n].s() = dimensionedScalar("0",dimMoles/dimVolume,Ci[n].C0membrane);
        forAll(cellAddressListFluid,celli)
        {
            Ci[n].s()[cellAddressListFluid[celli]] = Ci[n].C0aq;
        }
    } 

    Ci[n].s().write();
}
//**********************************************
//End loop
//*********************************************

//******************************************
//Set color function "alpha"
//alpha (membrane = 1)
//*****************************************
Info << "Reading Alpha\n";
volScalarField alpha(
    IOobject(
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimless, 1),
    zeroGradientFvPatchScalarField::typeName);

forAll(cellAddressListFluid,celli)
{
    alpha[cellAddressListFluid[celli]] = 0.0;
}


alpha.correctBoundaryConditions();


#include "createFieldsElectroStatic.H"