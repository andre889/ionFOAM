/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------


\*---------------------------------------------------------------------------*/
//#define Ampere 1

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
#include "../classes/timeController/electricControl.H"
#include "../classes/specie/specie.H"
#include "../classes/myTimeSchemes/myddtScheme.H"

#include "regionProperties.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define MYDEBUG 0
//#define MEMBRANE 1

#include <vector>
class material
{
private:
    /* data */
public:
    dimensionedScalar epsilon; //permittivity;
    dimensionedScalar mu; //viscosity;
    material(const dictionary & dict_);
    ~material();
};

material::material(const dictionary & dict_)
:   epsilon(dict_.lookup("er")),
    mu(dict_.lookup("mu"))
{
}

material::~material()
{
}


inline void updateQ(volScalarField& q, PtrList<specie>& Ci)
{       
                    q = q*0.0;
                    forAll(Ci,i)
                    {
                        q += Ci[i].s()*Ci[i].z*constant::physicoChemical::F;
                    }
                    q.correctBoundaryConditions();
};
 template<class TYPE>
    tmp<GeometricField<TYPE, fvPatchField, volMesh>> globalMeshToFluid(std::vector<label>& cellAddressListFluid,const GeometricField<TYPE, fvPatchField, volMesh> &vf,const fvMesh& mesh2)
    {
        dimensioned<TYPE> myZERO("ZERO",vf.dimensions(),TYPE());

        tmp<GeometricField<TYPE, fvPatchField, volMesh>> tnew(new GeometricField<TYPE, fvPatchField, volMesh>(IOobject(
                                                            "newf",
                                                            vf.time().timeName(),
                                                            mesh2,
                                                            IOobject::NO_READ,
                                                            IOobject::NO_WRITE),
                                                            mesh2,
                                                            myZERO,
                                                            zeroGradientFvPatchField<TYPE>::typeName));
        GeometricField<TYPE, fvPatchField, volMesh> & tX= tnew.ref();
        forAll(cellAddressListFluid,celli)
        {
            tX[celli] = vf[cellAddressListFluid[celli]];
        }

        tX.correctBoundaryConditions();
        return tnew;
    }
inline void updateFe(std::vector<label>& cellAddressListFluid, volVectorField& Fe, volVectorField& E, volScalarField& q, dimensionedScalar& rhof)
{
               forAll(cellAddressListFluid,celli)
                {
                        Fe[celli] = E[cellAddressListFluid[celli]]*q[cellAddressListFluid[celli]]/rhof.value(); //columbic force
                }
};

//Foam::debug = true;
int main(int argc, char *argv[])
{

std::vector<double> allDts;
allDts.reserve(10000);

//#include "postProcess.H"
#include "setRootCaseLists.H"
#include "createTime.H"  
#include "createFields.H"




#include "initContinuityErrs.H"

IOMRFZoneList MRF(meshFluid);

    Info << "....Note about non dimensionalizing the nernst plank poisson equation.... " << endl;
    Info << " for the equation: dn/dt = div(D div(n) + z*n*F* D/R/T div(phi)) and laplace(phi) = -z*n*F/epsilon" << endl;
    Info << " to the equation: dn/dt = div(div(n) + z*n*div(phi)) and laplace(phi) = -z*n" << endl;
    Info << " z* = RT/F, D* = 1, epsilon = 1/(RT)"  << endl;
    Info << "Diff " << Ci[0].Daq << endl;
    Info << " zF/epsilon " << Ci[0].z*constant::physicoChemical::F.value()/max(epsilon).value() << endl;
    Info << " z*F*D/R/T " << Ci[0].z * constant::physicoChemical::F.value() * Ci[0].Daq / constant::physicoChemical::R.value() / T.value() << endl;

    autoPtr<myddtScheme> timeSchemePtr(myddtScheme::New(mesh, mesh.schemesDict()));

    electricControl nppPISO(piso.dict(),mesh);

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    bool solveFluid = piso.dict().lookupOrDefault<bool>("solveFluid",true);

 

    //scalar newDT = runTime.deltaT().value();
    

    //Not sure if strictly needed, adding to ensure initialization of oldTime ptr within the GeometricField class before first loop.
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

            //Foam::Info << "what 1 " << Foam::tanh(Foam::mag(-3.0-runTime.value())) << Foam::endl;
            //scalar value = -2.0*(1.0-Foam::pow(Foam::tanh(Foam::mag(-3.0-runTime.value())),2));
            //Foam::Info << "what 1 " << value << Foam::endl;

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
                    // Info << "Step 1: Predict u*\n";
                    // Calc Fe(C_0,V_0)

                    updateQ(q, Ci);
                    E = -fvc::grad(voltage);
                    updateFe(cellAddressListFluid, Fe, E, q, rhof);


                    // Residual control for now
                     while(nppPISO.correct())//presidual > 1e-5)
                     {
                        #include "Ueqn.H" //move 08/30/2022

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
                           // UEqn = fvVectorMatrix(timeSchemePtr->ddt(U) + fvm::div(phi, U) +  turbulence->divDevSigma(U) + MRF.DDt(U));
                         }

                         updateQ(q, Ci);
                         E = -fvc::grad(voltage);
                         if (solveFluid)
                         {
                             updateFe(cellAddressListFluid, Fe, E, q, rhof);
                         }
                     }

                     Info << "errors\n";
                     {
                         volScalarField contErr(fvc::div(phi));

                         scalar sumLocalContErr = runTime.deltaTValue() *
                                                  mag(contErr)().weightedAverage(meshFluid.V()).value();

                         scalar globalContErr = runTime.deltaTValue() *
                                                contErr.weightedAverage(meshFluid.V()).value();
                         cumulativeContErr += globalContErr;

                         Info << "time?? step continuity errors : sum local = " << sumLocalContErr
                              << ", global = " << globalContErr
                              << ", cumulative = " << cumulativeContErr
                              << endl;
                     }
                    //Calculate and store error for all time varying functions...
                            Info << "last RK\n";

            } // END RK sub Loops
                                        Info << "after RK\n";

                UGlobal = fvc::reconstruct(phiGlobal); //added for debug only
                UGlobal.correctBoundaryConditions();  //added for debug only
                                        Info << "before write \n";

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

            allDts.push_back(runTime.deltaT().value());
        // runTime.setDeltaT(newDT);

        //Check mesh. If need refinement stop and write mesh
        

        } //END TIME


{
    OFstream osError("ci_end.csv", IOstream::ASCII, 2, IOstream::UNCOMPRESSED, false);
    osError << "x, c,\n";
    forAll(Ci[0].s(),celli)
    {
        scalar xloc = mesh.C()[celli].x();
        osError
            << xloc << "," << Ci[0].s()[celli] << ",\n";
    }
}

{
    OFstream osError("dt.csv", IOstream::ASCII, 2, IOstream::UNCOMPRESSED, false);
    osError << "i, dt,\n";
    forAll(allDts,i)
    {
        osError
            << i << "," << allDts[i] << ",\n";
    }

}

//Be Careful. if rp is greater than 1, then meshfluid is a new allocated object and needs deletion. 
if(rp.size() > 1)
{
    delete meshFluidPtr;
}


    Info<< "End\n" << endl;



    return 0;
}


// ************************************************************************* //
