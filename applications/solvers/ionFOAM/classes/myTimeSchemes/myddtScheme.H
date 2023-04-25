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

  virtual tmp<fvMatrix<scalar>>
   ddt
    (const GeometricField<scalar, fvPatchField, volMesh>& vf
    ) = 0;

 virtual tmp<fvMatrix<vector>>
   ddt
    (const GeometricField<vector, fvPatchField, volMesh>& vf
    ) = 0;
   
     virtual scalar getGamma()
     {
         return 1.0;
    }
    /// Access functions
    virtual label loopCount() = 0;
   
    /// Friend operators
    virtual label getImplicitStages() = 0;
    virtual scalar getAFunction(label i, label j) = 0;
    virtual bool loop() {
        currentLoop++;
        if (currentLoop == 2)
        {
            currentLoop = 0;
            // runTime.endSubCycle();
            return false;
        }
        return true;
    };
    virtual scalar getdTModifier(label i){return 1.0;};//modifier must be between 0 and 1 such that for loop i, t_i = t0 + dt*modifier;
    virtual scalar getBhatFunction(label i) { return 0; };
    virtual double getError() {return -1.234;};
};

}
#endif
