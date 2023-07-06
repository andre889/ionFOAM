        forAll(Ci,i)
        {
            volScalarField& ci = Ci[i].s();
            ci.storePrevIter();
            volScalarField Dfield(
                IOobject(
                "alpha",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
                mesh,
                dimensionedScalar("D",dimViscosity,Ci[i].Dmem),
                zeroGradientFvPatchScalarField::typeName); 
      

            surfaceScalarField Dsurface = diffusion ?
            fvc::interpolate(Dfield, mesh.interpolationScheme("interp(D)"))
            : 0.0*fvc::interpolate(Dfield, mesh.interpolationScheme("interp(D)"));

            surfaceScalarField theta = chemicalPotential ?
            (Dsurface * (1.0 - Ci[i].K)) / (fvc::interpolate(alpha) + (1.0 - fvc::interpolate(alpha)) * Ci[i].K) * fvc::snGrad(alpha) * mesh.magSf()
            : 0.0*(Dsurface * (1.0 - Ci[i].K)) / (fvc::interpolate(alpha) + (1.0 - fvc::interpolate(alpha)) * Ci[i].K) * fvc::snGrad(alpha) * mesh.magSf();

            surfaceScalarField eFlux = electroPhoretic ? 
            Ci[i].z * constant::physicoChemical::F* fvc::interpolate(Dfield, mesh.interpolationScheme("interp(D)")) * (linearInterpolate(E) & mesh.Sf())/constant::physicoChemical::R/T 
            : 0.0* constant::physicoChemical::F* fvc::interpolate(Dfield, mesh.interpolationScheme("interp(D)")) * (linearInterpolate(E) & mesh.Sf())/constant::physicoChemical::R/T; 
         

                fvScalarMatrix cEqn
                (
                    //fvm::ddt(Ci[i].s())
                    timeSchemePtr->ddt(ci)
                    - fvm::laplacian(Dsurface,ci) 
                    + fvm::div(theta,ci,"thetaScheme") 
                    + fvm::div(eFlux,ci) 
                    + fvm::div(phiGlobal,ci)
                );

                if(nppPISO.predictorElectric())
                {
                    nppPISO.storeErrorElectirc(cEqn.solve(mesh.solverDict("C")).initialResidual());
                }
                else
                {
                    scalar cEqnRelax = readScalar(piso.dict().lookup("cEqnRelax"));
                    cEqn.relax(cEqnRelax);
                    nppPISO.storeErrorElectirc(cEqn.solve(mesh.solverDict("Ccorr")).initialResidual());
                    scalar ciRelax = readScalar(piso.dict().lookup("ciRelax"));
                    Ci[i].s().relax(ciRelax);
                }

                    Ci[i].s().correctBoundaryConditions(); 
        }

    
