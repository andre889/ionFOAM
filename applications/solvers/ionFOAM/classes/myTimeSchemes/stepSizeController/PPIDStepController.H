#ifndef PPIDStepController_H
#define PPIDStepController_H

// dt^(n+1)=k(dt)^n(epsilon/(delta(n+1)))^alpha*(delta(n)/epsilon)^beta*(epsilon/delta(n-1))^gamma*(dt(n)/dt(n-1))^a*(dt(n-1)/dt(n-2))^b;
#include "fvCFD.H"

#include "dictionaryEntry.H"

#include "timeStepController.H"

namespace Foam
{

class PPIDStepController:
 public timeStepController
{

  private:
  
 scalar oldRelativeError = 1.0;
scalar oldOldRelativeError = 1.0;
 scalar oldoldDt = 1.0; 
 label counts = 0;
  public:
  
    // Declare Runtime type infomation
    TypeName("PPIDStepController");
    

        
    // Constructors
   PPIDStepController
    (
        const fvMesh& mesh,
        const dictionary& dict,
        scalar p
    );
    
   

    ~PPIDStepController();

   
    /// Access functions
            scalar newDt(scalar error){
              counts++;
            scalar dt = runTime_.deltaT().value(); //current dt value
            alpha = 6.0/(20.0*phat);
            scalar beta = -1.0/(20.0*phat);
            scalar a = 1.0;
            scalar gamma = -5/(20*phat);

              //return K*dt*Foam::pow(epsilon/deltaNp1,alpha)*Foam::pow(deltaN/epsilon,beta)*Foam::pow(epsilon/deltaNm1,gamma)*Foam::pow(dt/dtNm1,a)*Foam::pow(dtNm1/dtNm2,b);
             
            relativeError = error/epsilon;
            Info << "PPID Controller:\tEpsilon/error " << relativeError << endl;
 
            scalar dtToReturn =   K*dt*Foam::pow(1.0/relativeError,alpha)*Foam::pow(oldRelativeError,beta)*Foam::pow(1.0/oldOldRelativeError,gamma)*Foam::pow(dt/oldoldDt,a);
            
            if(counts < 3) //simple controller until 2 iterations have past. 
            {
               alpha = 1./(phat+1.0); //simple controller 
               dtToReturn =   K*dt*Foam::pow(1.0/relativeError,alpha);
            }

            Info << "PPID Controller:\tdt return: " << dtToReturn << endl;

            Info << "PPID Controller:\toldRelativeError: " << oldRelativeError << endl;
            Info << "PPID Controller:\toldOldRelativeError: " << oldOldRelativeError << endl;
            Info << "PPID Controller:\toldoldDt: " << oldoldDt << endl;
            
            oldoldDt = dt;
            oldOldRelativeError = oldRelativeError;
            oldRelativeError = relativeError;
            return dtToReturn; //*Foam::pow(deltaN/epsilon,beta)*Foam::pow(epsilon/deltaNm1,gamma)*Foam::pow(dt/dtNm1,a)*Foam::pow(dtNm1/dtNm2,b);
          };

};

}
#endif
