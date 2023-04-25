#ifndef timeStepController_H
#define timeStepController_H

#define LargeChangeInDt (oldDt-myNewDT)/oldDt > .1


// dt^(n+1)=k(dt)^n(epsilon/(delta(n+1)))^alpha*(delta(n)/epsilon)^beta*(epsilon/delta(n-1))^gamma*(dt(n)/dt(n-1))^a*(dt(n-1)/dt(n-2))^b;
#include "fvCFD.H"

#include "dictionaryEntry.H"

#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

class timeStepController
{
private:
    
  //- Disallow default bitwise copy construct
  timeStepController(const timeStepController&);
  //- Disallow default bitwise assignment
  void operator=(const timeStepController&);
  
  protected:
  
  //References
  const fvMesh& mesh_;
  const dictionary& dict_;
  Time &runTime_;

  public:
          scalar K = 0.95;
          scalar phat;
          scalar deltaNp1; //error for current step;
          scalar deltaN; // previous error
          scalar deltaNm1; // previous (x2) error
          scalar epsilon = 1e-5; //error tolerance
          scalar minDt = 1e-15; // minimum dt to avoid crashes. 
          scalar maxDt = 1e3;
          scalar alpha;
          scalar beta;
          scalar gamma;
          scalar a;
          scalar b;
          scalar dtNm1;
          scalar dtNm2;
          scalar relativeError; 
          scalar maxRelativeError = 10;
          scalar dtLimitFactor = 5;
    // Declare Runtime type infomation
    TypeName("default");
    

//declare table
    declareRunTimeSelectionTable(
        autoPtr,
        timeStepController,
        Normal,
        (
        const fvMesh& mesh,
        const dictionary& dict,
        scalar p
        ),
        (mesh,dict,p)
        );
        
    // Constructors
   timeStepController
    (
        const fvMesh& mesh,
        const dictionary& dict,
        scalar p
    );
    
   //reference to selected model
        static autoPtr<timeStepController> New
        (
              const fvMesh& mesh,
              const dictionary& dict,
              scalar p
        );
   
    
   //Destructor
   virtual ~timeStepController();

   
    /// Access functions
   //virtual double newDt();

        //Returns true if dt update was successful. 
        bool updateDT(scalar error)
          {
                scalar oldDt = runTime_.deltaT().value();
                scalar myNewDT =  this->newDt(error); 
                Info << "TimeStepController: new DT: " << myNewDT << "\n";
                //limit dt
                myNewDT = Foam::min(myNewDT, dtLimitFactor*oldDt);
                myNewDT = Foam::max(myNewDT, oldDt/dtLimitFactor);
                myNewDT = Foam::max(myNewDT,minDt);
                myNewDT = Foam::min(myNewDT,maxDt);
                //updateDt
                runTime_.setDeltaT(myNewDT);
                 Info << "TimeStepController: new DT Filtered: " << runTime_.deltaT().value() << "\n";
                bool HighRelativeError  = relativeError > maxRelativeError;
                bool stopConditions = (HighRelativeError && myNewDT != minDt);
                return !stopConditions;
          };

          virtual scalar newDt(scalar error){Info << "\n\n\n\nWARNING!!!\n\n\n\n";return runTime_.deltaT().value();};

};

}
#endif
