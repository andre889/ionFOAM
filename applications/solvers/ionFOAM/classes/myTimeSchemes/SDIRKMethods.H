// Diagonally Implicit Runge-Kutta methods for ordinary differential equations, a review CA Kennedy, MH Carpenter - 2016
// 5.1.1 ESDIRK
// Also see Diagonally implicit Runge–Kutta methods for stiff ODEs
//Christopher A.Kennedya,∗, Mark H.Carpenterb


//Time scheme controller 
// dt^(n+1)=k(dt)^n(epsilon/(delta(n+1)))^alpha*(delta(n)/epsilon)^beta*(epsilon/delta(n-1))^gamma*(dt(n)/dt(n-1))^a*(dt(n-1)/dt(n-2))^b;

#ifndef SDIRKMethods_H
#define SDIRKMethods_H

#define LoopStart currentLoop == 0
#define LoopEnd currentLoop == numStages + 1

#include "fvCFD.H"
#include "myddtScheme.H"
#include "scalar.H"
#include "scalarList.H"
#include "time.h"
#include "stepSizeController/timeStepController.H"

namespace Foam
{

class SDIRKMethods:
public myddtScheme
{
  private:
  Time &runTime_;
  scalar absoluteRelativeErrorSwitch;
  protected:

      label numStages; // s in paper reflects number of stages +1
      
      scalarListList a; // rk butcher table

      scalarList b;

      scalarList bHat;

      scalarList c; // step-fraction from 0 to 1

      scalar gamma;

      List<word> fieldList; //non unique names will break this program
      PtrList<PtrList<IOobject>> fList;

      scalar orignalTime;
      scalar orignalDt;
      label currentTimeLabel;
      
      autoPtr<timeStepController> myDtControllerPtr;

  public:
      // Constructor
      SDIRKMethods(
          const fvMesh &mesh,
          const dictionary &dict);

      // Destructor
      ~SDIRKMethods();

    
      tmp<fvMatrix<scalar>>
      ddt(const GeometricField<scalar, fvPatchField, volMesh> &vf);

         
      tmp<fvMatrix<vector>>
      ddt(const GeometricField<vector, fvPatchField, volMesh> &vf);

      template <class TYPE>
    void initializeRKMethod(const GeometricField<TYPE, fvPatchField, volMesh> &vf);

      //tmp<fvMatrix<vector>>
      //ddtv(const GeometricField<vector, fvPatchField, volMesh> &vf);
 template<class TYPE>
    void findFieldPtr(bool &found, label &foundIndex, const GeometricField<TYPE, fvPatchField, volMesh> &vf);
 
 template<class TYPE>
    tmp<GeometricField<TYPE, fvPatchField, volMesh>> getX(const GeometricField<TYPE, fvPatchField, volMesh> &vf,label index);

 template<class TYPE>
    void updateFFromPreviousCalc(label vfIndex ,const GeometricField<TYPE, fvPatchField, volMesh> &vf);

template<class TYPE>
    double findError(GeometricField<TYPE, fvPatchField, volMesh>& vf); //fieldList[fieldi],fList[fieldi]

    //Add field to list so that quantities such as f can be allocated.
template<class TYPE>
    void allocateF(label& foundIndex, const GeometricField<TYPE, fvPatchField, volMesh> &vf);

      template<class TYPE>
      void loopRestart();
      
      label loopCount()
      {
          return numStages;
      };
      
      void loopRestart()
      {
          currentLoop = 0; 
         forAll(fieldList,fieldi)
         {
            Foam::HashTable<Foam::regIOobject*>::const_iterator fieldiPtr = fList[fieldi][0].db().find(fieldList[fieldi]);
                          IOobject& object1 = static_cast<IOobject&>(*fieldiPtr());
                        //Info << "which is a   " << object1.type() << "\n";
                        //Info << "which is a   " << object1.headerClassName() << "\n";
                        
                      
                    if(object1.type() == "volScalarField")
                    {
                              volScalarField& vfFieldi = static_cast<volScalarField&>(object1);
                              vfFieldi = vfFieldi.oldTime();
                         }
                         else if(object1.type() == "volVectorField")
                         {
                              volVectorField& vfFieldi = static_cast<volVectorField&>(object1);
                              vfFieldi = vfFieldi.oldTime();
                         }
                         else
                         {
                        
                            Info << "\n\nSorry but type " << object1.type()  << " is not supported!\n";
                            
                        // code block
                    }
          
         }
         clearFPtr();

         //Reset time as if it had progressed to t0+newDeltat;
        scalar theNewDt = runTime_.deltaT().value();
        scalar t0 = (orignalTime - orignalDt);
        runTime_.setTime(t0+theNewDt, currentTimeLabel);

      };

      bool loop()
      {

          if(LoopStart)
          {
              storeCurrentTimeValues();
          }
         
          currentLoop++;
         
          if (LoopEnd)
          {

              if(myDtControllerPtr.valid()) // if ptr exists means there is a Time controller to use
              {
                scalar error = gatherError(); 

                 Info << "SDIRK:\tThe Error estimate is " << error << "\n"; // Debug
                
                
                if(myDtControllerPtr->updateDT(error))
                {
                    Info << "SDIRK:\tSuccessful new dt = " << runTime_.deltaT().value() << endl;
                }
                else
                {
                    Info << "SDIRK:\tRestarting RK loop with dt = " <<  runTime_.deltaT().value() << endl;
                    loopRestart();
                    return loop();
                }
             }
             
              currentLoop = 0;
              clearFPtr();

              //Reset time

              runTime_.setTime(orignalTime, currentTimeLabel);

              return false;
          }
          scalar timeN = (orignalTime - orignalDt); // This is the time just before the implicit stage n+1
          scalar rkTime = timeN + getdTModifier(currentLoop - 1) * orignalDt;

          runTime_.setTime(rkTime, currentTimeLabel);
                                 
          return true;
      };

    scalar getGamma()
    {
        return gamma;
    };
    scalar getdTModifier(label i) { return c[i]; };
    scalar getAFunction(label i, label j) {
        //Info << "A is returned" << a[i][j] <<"\n";
        return a[i][j];
    };

    void clearFPtr()
    {
        forAll(fList, vfi){
            fList[vfi].clear();
            fList[vfi].setSize(numStages - 1);
        } 
    };

    //Stores the current Time value, dt and time label
    void storeCurrentTimeValues()
    {
              //Store Old time
              orignalTime = runTime_.value();
              //Store Old dt
              orignalDt = runTime_.deltaT().value();
              currentTimeLabel = runTime_.timeIndex();
    };

   
    label getImplicitStages() { return 2; }; // Number of "Rs" to store
    scalar getBhatFunction(label i) { return bHat[i]; };

   

    //For embeded methods, calculates the largest error (delta) from all variables in the system.
    scalar gatherError()
    {
        scalar error = VSMALL;
        forAll(fieldList,fieldi)
        {

            Foam::HashTable<Foam::regIOobject*>::const_iterator fieldiPtr = fList[fieldi][0].db().find(fieldList[fieldi]);
            IOobject& object1 = static_cast<IOobject&>(*fieldiPtr());

                if(object1.type() == "volScalarField")
                {
                        volScalarField& vfFieldi = static_cast<volScalarField&>(object1);
                        error = Foam::max(error,findError(vfFieldi));
                }
                else if(object1.type() == "volVectorField")
                {
                    volVectorField& vfFieldi = static_cast<volVectorField&>(object1);
                    error = Foam::max(error,findError(vfFieldi));
                }
                else
                {
                    Info << "\n\nSorry but type " << object1.type()  << " is not supported!\n";
                }
                
        
        }

        return error;
    };

};

}

#endif
