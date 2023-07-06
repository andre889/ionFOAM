//***********************************************
//  Date:  05/05/2023
//  Revision: 1.0
//  Author: Austin J. Andrews
//  Email: andre889@umn.edu
//  Affiliation: University of Minnesota Twin Cities, Department of Mechanical Engineering, Hogan Lab

//  Description: This class is a parent class for all of the new time scheme methods. This class mimics the use of the fvm::ddt (matrix construct) functions, but written in a way to accept multi-step time schemes.
//  This can potentially be merged into the existing usage of fvm::ddt, but has yet to be done.
// **********************************************

#ifndef myddtScheme_H
#define myddtScheme_H

#include "fvCFD.H"

#include "dictionaryEntry.H"

#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

class myddtScheme
{

  private:
  
  //- Disallow default bitwise copy construct
  myddtScheme(const myddtScheme&);
  //- Disallow default bitwise assignment
  void operator=(const myddtScheme&);
  
  protected:
  
  //References
  const fvMesh& mesh_;
  const dictionary& dict_;
  label currentLoop = 0;
  //Stored data
  public:
  
    // Declare Runtime type infomation
    TypeName("myddtScheme");
    

    //declare table
    declareRunTimeSelectionTable(
        autoPtr,
        myddtScheme,
        Normal,
        (
        const fvMesh& mesh,
        const dictionary& dict
        ),
        (mesh,dict)
        );
        
    // Constructors
   myddtScheme
    (
        const fvMesh& mesh,
        const dictionary& dict
    );
    
   //reference to selected model
        static autoPtr<myddtScheme> New
        (
              const fvMesh& mesh,
              const dictionary& dict
        );
   
    
   //Destructor
   virtual ~myddtScheme();

   //Templates disallow virtual functions, but this should be virtual.
   template <class TYPE>
   tmp<fvMatrix<TYPE>>
     ddt
    (
       const GeometricField<TYPE, fvPatchField, volMesh>& vf
    ){tmp<fvMatrix<TYPE>> tfvm(new fvMatrix<TYPE>(vf,vf.dimensions() * dimVol / dimTime)); return tfvm;};
   
    virtual scalar getGamma()
    {
         return 1.0;
    }

    /// Access functions
     label loopCount(){return 0;};
     label getImplicitStages(){return 0;};
     scalar getAFunction(label i, label j){return 0;};
     scalar getdTModifier(label i){return 1.0;};//modifier must be between 0 and 1 such that for loop i, t_i = t0 + dt*modifier;
     scalar getBhatFunction(label i) { return 0; };
     double getError() {return -1.234;};

     bool loop() {
        currentLoop++;
        if (currentLoop == 2)
        {
            currentLoop = 0;
            return false;
        }
        return true;
    };

};

}
#endif
