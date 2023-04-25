Info << "Create mesh for time = "
    << runTime.timeName() << nl << endl;

fvMesh mesh(
    IOobject(
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ));


Info << "Done creating mesh!\n";
Info << "Default Mesh has: " << mesh.C().size() << " cells\n";

fvMesh meshFluid(
    IOobject(
        "fluid",
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ));

Info << "Fluid Mesh has: " << meshFluid.C().size() << " cells\n";

Info << "Creating addressing for fluid to global mesh regions...\n";

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

    labelUList cellAddressListFluid = cellAddressing;
    labelUList boundaryRegionAddressList = boundaryRegionAddressing;

    labelList faceAddressList(faceAddressing1.size());
    for (label i = 0; i < faceAddressing1.size(); ++i)
    {
        faceAddressList[i] = faceAddressing1[i];
    }
      
    /*List<labelList> faceAddressListList;
    for (label i = 0; i != Pstream::nProcs(); ++i)
    {
        faceAddressListList.append(faceAddressList);
        //faceAddressListList.insert(i, new labelUList(faceAddressList));
    }

        labelUList boundaryRegionAddressList = boundaryRegionAddressing;
      
            for (label i = 0; i != faceAddressListList.size(); ++i)
            {
                Info << "Hello from processor " << Pstream::myProcNo() << "! I am working on "
                     << faceAddressListList[i] << endl;
            }
        

            Info << "Fluid to global mesh mapping complete. " << endl;

            Info << "Creating parallel mappings " << endl;

            labelListList faceProcAddressingList(Pstream::nProcs());

            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                Time database(
                    Time::controlDictName,
                    runTime.rootPath(),
                    runTime.caseName() / ("processor" + Foam::name(proci)));

                fvMesh meshToAdd(
                    IOobject(
                        "fluid",
                        database.timeName(),
                        database));
                labelIOList faceAddressingProci(
                    IOobject(
                        "polyMesh/faceProcAddressing",
                        database.time().constant(),
                        meshToAdd,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE));

                faceProcAddressingList[proci].resize(faceAddressingProci.size());
                for (label j = 0; j != faceProcAddressingList[proci].size(); ++j)
                {
                    faceProcAddressingList[proci][j] = faceAddressingProci[j];
                }

                Info << "Proc list " << proci << " \n"
                     << faceProcAddressingList[proci] << endl;
            }

            //faceProcAddressing
            
            */
            //exit(2);