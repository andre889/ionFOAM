/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
//***********************************************************
//  Date:  05/05/2023
//  Revision: 1.0
//  Author: Austin J. Andrews
//  Email: andre889@umn.edu
//  Affiliation: University of Minnesota Twin Cities, Department of Mechanical Engineering, Hogan Lab

//  Description: Based on the pisoFOAM solver, this program, coined "ionFoam", solves 
//  the NS-NPP set of equations with application to ion selective membranes (ISE) with
//  non-equal partitioning of ions at phase boundaries. More information can be found 
//  in the study titled "Determination of Forced Convection Effects on the Response of 
//  Ion Selective Electrodes via Numerical Solution to the Navier-Stokes-Nernst-Plank-Poisson 
//  Equations"
// **********************************************
\*---------------------------------------------------------------------------*/

#include <vector>

#include "fvCFD.H"
#include "viscosityModel.H"
#include "incompressibleMomentumTransportModels.H"
#include "pisoControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "IOobjectList.H"
#include "lduPrimitiveMesh.H"
#include "simpleMatrix.H"
#include "OFstream.H"
#include "regionProperties.H"

#include "../classes/timeController/electricControl.H"
#include "../classes/specie/specie.H"
#include "../classes/myTimeSchemes/myddtScheme.H"

#include "material.H"
#include "auxFunctions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#include "setRootCaseLists.H"
#include "createTime.H"  
#include "createFields.H"
#include "initContinuityErrs.H"

    IOMRFZoneList MRF(meshFluid);

    autoPtr<myddtScheme> timeSchemePtr(myddtScheme::New(mesh, mesh.schemesDict()));

    electricControl nppPISO(piso.dict(),mesh);

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    bool solveFluid = piso.dict().lookupOrDefault<bool>("solveFluid",true);
    

    //Ensure initialization of oldTime ptr within the GeometricField class before first loop.
    forAll(Ci,i)
    {
        Ci[i].s().oldTime();
    }

    U.oldTime();
    U.write();

    while (runTime.run())
    {

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        dimensionedScalar rDeltaTGamma("rDeltaTGamma", dimless / dimTime, 1.0 / (runTime.deltaT().value() * timeSchemePtr->getGamma()));
      
        //Store Old Times
        U.storeOldTime();

        forAll(Ci,i)
        {
            Ci[i].s().storeOldTime();
        }
        
       
        // Pressure-velocity PISO corrector
        // see https://doi.org/10.1016/0021-9991(86)90099-9 
        // and https://doi.org/10.1080/104077901753306601
            U.correctBoundaryConditions();
            p.correctBoundaryConditions();
            voltage.correctBoundaryConditions();
            forAll(Ci,i)
            {
                Ci[i].s().correctBoundaryConditions();
            }

            while (timeSchemePtr->loop())
            {

                    updateQ(q, Ci);
                    E = -fvc::grad(voltage);
                    updateFe(cellAddressListFluid, Fe, E, q, rhof);

                     while(nppPISO.correct())
                     {
                        #include "Ueqn.H" 

                        //store previous iterations for relaxation.
                        forAll(Ci,i)
                        {
                            Ci[i].s().storePrevIter();
                        }
                        voltage.storePrevIter();

                        #include "pEqn.H"

                         while (nppPISO.loopElectric())
                         {
                             voltage.storePrevIter();
                             updateQ(q, Ci);
#include "voltageEqn.H"
#include "ionEqn.H"
                         }
                         updateQ(q, Ci);
                         E = -fvc::grad(voltage);
                         if (solveFluid)
                         {
                             updateFe(cellAddressListFluid, Fe, E, q, rhof);
                         }
#include "pEqn.H" //Correct BodyForce Term
                         if (solveFluid && multiRegion)
                         {
#include "localToGlobalMapping.H"
                         }
#include "ionEqn.H"
                         if (solveFluid)
                         {
                            viscosity->correct();
                            turbulence->correct();
                         }

                         updateQ(q, Ci);
                         E = -fvc::grad(voltage);
                         if (solveFluid)
                         {
                             updateFe(cellAddressListFluid, Fe, E, q, rhof);
                         }
                     }

                     {
                         volScalarField contErr(fvc::div(phi));

                         scalar sumLocalContErr = runTime.deltaTValue() *
                                                  mag(contErr)().weightedAverage(meshFluid.V()).value();

                         scalar globalContErr = runTime.deltaTValue() *
                                                contErr.weightedAverage(meshFluid.V()).value();
                         cumulativeContErr += globalContErr;

                         Info << "time step continuity errors : sum local = " << sumLocalContErr
                              << ", global = " << globalContErr
                              << ", cumulative = " << cumulativeContErr
                              << endl;
                     }
                    Info << "End RK loop\n";
             } // END RK sub Loops
                                        

                UGlobal = fvc::reconstruct(phiGlobal); //added for debug only
                UGlobal.correctBoundaryConditions();  //added for debug only

            if(runTime.write())
            {
                //force write fluid data. Current bug with just fluid simulations.
                U.write();
                p.write();
                //turbuluence->write();
            }
        

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        } //END TIME Loop


//Be Careful. if rp is greater than 1, then meshfluid is a "new" allocated object and needs deletion. 
if(rp.size() > 1)
{
    delete meshFluidPtr;
}


    Info<< "End\n" << endl;



    return 0;
}


// ************************************************************************* //
