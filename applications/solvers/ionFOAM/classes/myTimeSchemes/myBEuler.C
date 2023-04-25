#include "myBEuler.H"
namespace Foam
{
    //Add to Selection Table
    defineTypeNameAndDebug(myBEuler, 0);
    addToRunTimeSelectionTable(myddtScheme, myBEuler, Normal);
    
    
myBEuler::myBEuler
    (
        const fvMesh& mesh,
        const dictionary& dict
    )
   :
    myddtScheme(mesh,dict)
    {
    
    }
  

 myBEuler::~myBEuler(){}


  tmp<fvMatrix<scalar>>
   myBEuler::ddt
    (const GeometricField<scalar, fvPatchField, volMesh>& vf
    )
    {
        dimensionedScalar rDeltaT("rDeltaT", dimless / dimTime, 1.0 / (vf.time().deltaT().value()));

        tmp<fvMatrix<scalar>> tfvm( new fvMatrix<scalar> (vf,
            vf.dimensions()*dimVol/dimTime));

    const fvMesh& mesh_ = vf.mesh();

    //For moving mesh Vsc() otherwise V()
    tfvm.ref().diag() = rDeltaT*mesh_.V();
    tfvm.ref().source() = rDeltaT*mesh_.V()*vf.oldTime().primitiveField();

    return tfvm;
    };

    tmp<fvMatrix<vector>>
   myBEuler::ddt
    ( const GeometricField<vector, fvPatchField, volMesh>& vf
    )
    {
        dimensionedScalar rDeltaT("rDeltaT", dimless / dimTime, 1.0 / (vf.time().deltaT().value()));
        tmp<fvMatrix<vector>> tfvm( new fvMatrix<vector> (vf,
            vf.dimensions()*dimVol/dimTime));

    const fvMesh& mesh_ = vf.mesh();

    //For moving mesh Vsc() otherwise V()
    tfvm.ref().diag() = rDeltaT*mesh_.V();
    tfvm.ref().source() = rDeltaT*mesh_.V()*vf.oldTime().primitiveField();

    return tfvm;
    };

}
