 

if(nppPISO.predictorElectric())
{
        Info << "Voltage Predictor Loop\n";
       
       
        volScalarField sumAn = 0.0*epsilon; 
        volScalarField sumBn = 0.0*q;

        dimensionedScalar effectiveDt = runTime.deltaT()*timeSchemePtr->getGamma();
        
        forAll(Ci, i)
        {
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

            sumAn += Ci[i].z*Ci[i].z* effectiveDt * constant::physicoChemical::F * constant::physicoChemical::F * (Dfield)*Ci[i].s() / constant::physicoChemical::R / T;

            surfaceScalarField Dsurface = diffusion ?
            fvc::interpolate(Dfield, mesh.interpolationScheme("interp(D)"))
            : 0.0*fvc::interpolate(Dfield, mesh.interpolationScheme("interp(D)"));

            surfaceScalarField theta = chemicalPotential ?
            (Dsurface * (1.0 - Ci[i].K)) / (fvc::interpolate(alpha) + (1.0 - fvc::interpolate(alpha)) * Ci[i].K) * fvc::snGrad(alpha) * mesh.magSf()
            : 0.0*(Dsurface * (1.0 - Ci[i].K)) / (fvc::interpolate(alpha) + (1.0 - fvc::interpolate(alpha)) * Ci[i].K) * fvc::snGrad(alpha) * mesh.magSf();
            
            volScalarField Hc =  fvc::div(theta,Ci[i].s(),"predictorMTScheme") - fvc::laplacian(Dsurface,Ci[i].s()); 

            sumBn += Ci[i].z * constant::physicoChemical::F * effectiveDt * Hc;
        }

  label nccEND = readLabel(piso.dict().lookup("nNonOrthogonalCorrectorsV"));
        volScalarField sourceTot = -q + sumBn;

        //no- orthogonal corrector loop
        for (label ncc = 0; ncc < nccEND; ncc++) 
        {   
            fvScalarMatrix vEqn
            (
                fvm::laplacian(epsilon + sumAn,voltage)  == sourceTot 
            );
        
            vEqn.setReference(voltageRefCell, voltageRefValue);
            if(voltage.needReference())
            {
                vEqn.source() = vEqn.source()  - gAverage(vEqn.source());
            }
            scalar error = vEqn.solve().initialResidual(); //Do not store the error for the "Predictor" equation. Error will be managed on the next iteration.
        }
}
else
{
        label nccEND = readLabel(piso.dict().lookup("nNonOrthogonalCorrectorsV"));

         volScalarField sourceTot = -q;

        for (label ncc = 0; ncc < nccEND; ncc++) 
        {
                fvScalarMatrix vEqn(
                    fvm::laplacian(epsilon, voltage) == -q);
                vEqn.setReference(voltageRefCell, voltageRefValue);
                if(voltage.needReference())
                {
                    vEqn.source() = vEqn.source()  - gAverage(vEqn.source());
                }
                scalar error = vEqn.solve().initialResidual();
                nppPISO.storeErrorElectirc(error);
        }
    voltage.relax();
}

voltage.correctBoundaryConditions();
//Update E
E = -fvc::grad(voltage);

