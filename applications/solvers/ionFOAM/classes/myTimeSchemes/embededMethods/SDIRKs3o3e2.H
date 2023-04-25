// https://doi.org/10.1093/imamat/24.3.293
// alexander's method https://www.jstor.org/stable/2156678

#ifndef SDIRKs3o3e2_H
#define SDIRKs3o3e2_H

#include "fvCFD.H"
#include "../myddtScheme.H"
#include "scalar.H"
#include "scalarList.H"
#include "../SDIRKMethods.H"
#include "../stepSizeController/timeStepController.H"
namespace Foam
{

class SDIRKs3o3e2:
public SDIRKMethods
{

  protected:

    
  private:

  public:

      // Define name of distribution
      TypeName("SDIRKs3o3e2");

      // Constructor
      SDIRKs3o3e2(
          const fvMesh &mesh,
          const dictionary &dict);

      // Destructor
      ~SDIRKs3o3e2();

   
        
      
};

}

#endif
