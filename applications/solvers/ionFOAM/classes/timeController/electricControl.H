//  Date:  05/05/2023
//  Revision: 1.0
//  Author: Austin J. Andrews
//  Email: andre889@umn.edu
//  Affiliation: University of Minnesota Twin Cities, Department of Mechanical Engineering, Hogan Lab

//  Description: This class manages the PISO algorithm as described in the Determination of Forced Convection Effects on the Response of Ion Selective Electrodes via Numerical Solution to the Navier-Stokes-Nernst-Plank-Poisson Equations
// **********************************************


class electricControl
{
private:
    bool electrostaticPredictor;
    label electricIter; //Count iterations in the NPP equation loop
    scalar maxErrorE; //Here error is defined as the l2norm of the residual, E is for the NPP equations.
    scalar maxError;
    const dictionary &dict_;
    const fvMesh &mesh_;

public:
    label outerPISOiter; // PISO corrector number
    electricControl(const dictionary &, const fvMesh &);
    ~electricControl();

    bool correct()
    {
        outerPISOiter++;

        if (outerPISOiter > 1)
        {
            Info << "PISO:\t\tmax outer loop error " << maxError << endl;
            electrostaticPredictor = false;
        }

        label nCorrectorsMax = dict_.lookupOrDefault<label>("nCorrectorsMax", 4);
        label nCorrectorsMin = dict_.lookupOrDefault<label>("nCorrectorsMin", 2);
        scalar errorTol = dict_.lookupOrDefault<scalar>("errorTol", 1e-6);

        // Limit minimum counts
        nCorrectorsMax = max(nCorrectorsMax, nCorrectorsMin);

        bool minCountsReached = outerPISOiter >= nCorrectorsMin;
        bool maxCountsExceeded = outerPISOiter >= nCorrectorsMax;
        bool errorWithinTolerance = (maxError < errorTol);

        bool loopStop = (errorWithinTolerance && minCountsReached) || maxCountsExceeded;

        bool loopContinue = !loopStop; // return if loop should continue
        if (!errorWithinTolerance && maxCountsExceeded)
        {
            Info << endl
                 << "DANGER! Error not meet within " << nCorrectorsMax << " iterations!"
                 << " Consider increasing max iterations." << endl
                 << endl;
        }

        if (outerPISOiter == 1)
        {
            loopContinue = true; // Force First iter to always run loop.
        }
        // Reset error.

        maxError = VSMALL;

        if (!loopContinue)
        {
            reset();
        }
        else
        {
            Info << "PISO:\t\t outer loop iteration:  " << outerPISOiter << endl;
            resetElectric();
        }

        return loopContinue;
    };

    bool loopElectric()
    {

        electricIter++;

        label nNPPCorrectorsMax = dict_.lookupOrDefault<label>("maxNPPCorrectors", 4);
        label nNPPCorrectorsMin = dict_.lookupOrDefault<label>("minNPPCorrectors", 2);
        scalar errorTol = dict_.lookupOrDefault<scalar>("errorTol", 1e-6);

        // Limit minimum counts
        nNPPCorrectorsMax = max(nNPPCorrectorsMax, nNPPCorrectorsMin);

        bool minCountsReached = electricIter > nNPPCorrectorsMin;
        bool maxCountsExceeded = electricIter >= nNPPCorrectorsMax;
        bool errorWithinTolerance = (maxErrorE < errorTol);

        bool loopStop = (errorWithinTolerance && minCountsReached) || maxCountsExceeded;
        bool loopContinue = !loopStop; 

        if (!errorWithinTolerance && maxCountsExceeded)
        {
            Info << endl
                 << "DANGER! INNER LOOP: Error not meet within " << nNPPCorrectorsMax << " iterations!"
                 << " Consider increasing max iterations." << endl
                 << endl;
        }

        if (electricIter == 1)
        {
            loopContinue = true; 
        }

        if (!loopContinue)
        {
            Info << "PISO:\t\t\t****final inner loop error " << maxErrorE << endl;
        }
        else
        {
            if (electricIter > 1)
            {
                Info << "PISO:\t\t\tmax inner loop error " << maxErrorE << endl;
            }
            Info << "PISO:\t\t\tinner loop " << electricIter << "\n";
        }

        maxErrorE = VSMALL; 

        return loopContinue;
    };

    void storeErrorElectirc(scalar error)
    {
        maxErrorE = max(error, maxErrorE);
        maxError = max(error, maxError);
    };

    void storeError(scalar error)
    {
        maxError = max(error, maxError);
    };

    bool predictorElectric()
    {
        bool loopContinue = electrostaticPredictor && electricIter == 1 && outerPISOiter == 1;
        return loopContinue;
    };

    void resetElectric()
    {
        electricIter = 0;
        electrostaticPredictor = dict_.lookupOrDefault<bool>("NPPPredictor", true);
        maxErrorE = 1e9;
    };
    void reset()
    {
        resetElectric();
        outerPISOiter = 0;
        maxError = VSMALL;
    };
};

electricControl::electricControl(const dictionary &dict, const fvMesh &mesh)
    : dict_(dict),
      mesh_(mesh)
{
    reset();
}

electricControl::~electricControl()
{
}
