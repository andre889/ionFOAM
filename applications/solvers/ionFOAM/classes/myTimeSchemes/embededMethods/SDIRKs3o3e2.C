#include "SDIRKs3o3e2.H"
namespace Foam
{
    //Add to Selection Table
    defineTypeNameAndDebug(SDIRKs3o3e2, 0);
    addToRunTimeSelectionTable(myddtScheme, SDIRKs3o3e2, Normal);
    
    //x≈0.43586652150845899941601945119355684
SDIRKs3o3e2::SDIRKs3o3e2
    (
        const fvMesh& mesh,
        const dictionary& dict
    )
   :
    SDIRKMethods(mesh,dict)
    {
         scalar p = 3; // it is a 3rd order method.
         bool useController = dict.lookupOrDefault<bool>("variableTime",true);
         if(useController)
         {
              myDtControllerPtr = autoPtr<timeStepController>(timeStepController::New(mesh, dict,p));
         }

        numStages = p; // s in paper reflects number of stages +1
        gamma = 0.43586652150845899941601945119355684; //root of  x^3 - 3 x^2 + 3×x/2 - 1/6 = 0
        scalar alpha = gamma;
        scalar tau = (alpha*alpha - 3.0*alpha/2.0 + 1.0/3.0)/(alpha*alpha - 2.0*alpha + 0.5);
        scalar b1 = (tau/2.0 - 1.0/6.0)/((tau-alpha)*(1.0-alpha));
        scalar b2 = (alpha/2.0 - 1.0/6.0)/((alpha-tau)*(1.0-tau));
        scalar c1 = (tau - 0.5)/(tau - alpha);
        scalar c2 = (alpha-0.5)/(alpha-tau);

        b.setSize(numStages);
        
        b[0] = (1.0-gamma);
        b[1] = gamma;
        b[2] = 0;

        a.setSize(numStages);
        forAll(a, i)
        {
            a[i].setSize(numStages);
        }
     
        a[0][0] = gamma;
        a[0][1] = 0;
        a[0][2] = 0;

        a[1][0] = (tau-gamma);
        a[1][1] = gamma;
        a[1][2] = 0;

        a[2][0] = b1;
        a[2][1] = b2;
        a[2][2] = gamma;

        bHat.setSize(2);


        bHat[0] = c1;
        bHat[1] = c2;
        bHat[2] = 0;

        c.setSize(3);


        c[0] = gamma;
        c[1] = tau;
        c[2] = 1.0;
    };

    SDIRKs3o3e2::~SDIRKs3o3e2() {  };


}
