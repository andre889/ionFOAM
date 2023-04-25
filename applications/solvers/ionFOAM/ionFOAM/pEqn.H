if (solveFluid)
{
        //Solve P*(U*,C,V,u(for Fe))
        //fvVectorMatrix UEqn = UEqnNoTurb + turbulence->divDevReff(U) + MRF.DDt(U) ;
        scalar currentError = 1e31;
        scalar errorTolerance = piso.dict().lookupOrDefault<scalar>("fluidTolerance",1e-6);
        label nFluidLoops = piso.dict().lookupOrDefault<label>("nFluidLoops",1);

        label numCounts = 0;
        while(currentError > errorTolerance && numCounts < nFluidLoops)
        {
        numCounts++;

        volScalarField rAU(1.0 / UEqn.A());
        volVectorField HbyA(constrainHbyA(rAU * (UEqn.H() + Fe), U, p)); //Makes HbyA consistant with BCs.... should include FebyA?
        //+ fvc::interpolate(rAU)*fvc::ddtCorr(U, phi) //is this RHS? Not needeD? https://www.sciencedirect.com/science/article/abs/pii/S0045793014000334?via%3Dihub 
        surfaceScalarField phiHbyA(
            "phiHbyA",
            fvc::flux(HbyA) //+ fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
        );
        MRF.makeRelative(phiHbyA);

        adjustPhi(phiHbyA, U, p); //hmm some kind of incompressible correction?
        //Add RHS contributions

        //surfaceScalarField phiHbyAPlusFe = FebyA+phiHbyA;

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAU, MRF);
        //\nabla P = (\phi_hbya - \phi)a_s/Area;

        // Non-orthogonal pressure corrector loop
        label nccEND = readLabel(piso.dict().lookup("nNonOrthogonalCorrectors"));
        for (label ncc = 0; ncc < nccEND; ncc++) //while (piso.correctNonOrthogonal())
        {
            // Pressure corrector

            fvScalarMatrix pEqn(
                fvm::laplacian(rAU, p) == fvc::div(phiHbyA));

            pEqn.setReference(pRefCell, pRefValue); //p.select(piso.finalInnerIter())
            scalar error = pEqn.solve(mesh.solverDict("p")).initialResidual();
            nppPISO.storeError(error);

            if (nccEND - 1 == ncc)
            {
                Info << "phi correction\n";
                currentError = error;
                //By haveing phiHbyA on face Rhie and Chow correction is satisifed!
                phi = phiHbyA - pEqn.flux();
            }
        }

        //p.relax();
        
        //U = (rAU+rDIAGUT)*(A*phi - fvc::grad(p) + sourceDUDT); //rhie chow time correction?
        U = HbyA - rAU * fvc::grad(p); //Dont think the extra term is needed as Fe is in source term in Ueqn.
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();
        }
}
