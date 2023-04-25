
#ifndef myBEuler_H
#define myBEuler_H

#include "fvCFD.H"
#include "myddtScheme.H"
namespace Foam
{

class myBEuler:
public myddtScheme
{

  protected:

    
  private:
  
  public:
  
   //Define name of distribution
   TypeName("myBEuler");
   
   // Constructor
   myBEuler
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

   
   //Destructor
   ~myBEuler();
   
    tmp<fvMatrix<scalar>>
   ddt
    (const GeometricField<scalar, fvPatchField, volMesh>& vf
    );

    tmp<fvMatrix<vector>>
   ddt
    (const GeometricField<vector, fvPatchField, volMesh>& vf
    );
    
    label loopCount()
    {
        return 1;
    };
    label getImplicitStages() { return 0; };
    scalar getAFunction(label i, label j) { return 0; };
};

}

#endif
