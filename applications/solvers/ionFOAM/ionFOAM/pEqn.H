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
        volVectorField HbyA(constrainHbyA(rAU * (UEqn.H() + Fe), U, p)); //+ fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)  not needed https://www.sciencedirect.com/science/article/abs/pii/S0045793014000334?via%3Dihub 
        surfaceScalarField phiHbyA(
            "phiHbyA",
            fvc::flux(HbyA) //+ fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
        );
        MRF.makeRelative(phiHbyA);

        adjustPhi(phiHbyA, U, p);
       

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAU, MRF);

        // Non-orthogonal pressure corrector loop
        label nccEND = readLabel(piso.dict().lookup("nNonOrthogonalCorrectors"));
        for (label ncc = 0; ncc < nccEND; ncc++) 
        {
            // Pressure corrector
            fvScalarMatrix pEqn(
                fvm::laplacian(rAU, p) == fvc::div(phiHbyA));

            pEqn.setReference(pRefCell, pRefValue); 
            scalar error = pEqn.solve(mesh.solverDict("p")).initialResidual();
            nppPISO.storeError(error);

            if (nccEND - 1 == ncc)
            {
                Info << "phi correction\n";
                currentError = error;
                //By having phiHbyA on faces Rhie and Chow correction is satisifed. This is consistant with the general PISO algorithm provided by OpenFOAM's pisoFOAM.
                phi = phiHbyA - pEqn.flux();
            }
        }

            U = HbyA - rAU * fvc::grad(p);
            U.correctBoundaryConditions();
            p.correctBoundaryConditions();
        }
}
