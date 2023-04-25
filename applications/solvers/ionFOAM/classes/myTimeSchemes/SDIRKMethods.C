#include "SDIRKMethods.H"
namespace Foam
{

    SDIRKMethods::SDIRKMethods(
        const fvMesh &mesh,
        const dictionary &dict)
        :
        runTime_(const_cast<Time &>(mesh.time())),
        myddtScheme(mesh, dict)
        {
            absoluteRelativeErrorSwitch = dict.lookupOrDefault<scalar>("absoluteRelativeErrorSwitch",1);

        };

    SDIRKMethods::~SDIRKMethods(){};

    tmp<fvMatrix<scalar>>
    SDIRKMethods::ddt(const GeometricField<scalar, fvPatchField, volMesh> &vf)
    {
        dimensionedScalar rDeltaTGamma("rDeltaTGamma", dimless / dimTime, 1.0 / (vf.time().deltaT().value() * this->getGamma()));

        initializeRKMethod(vf);

        tmp<fvMatrix<scalar>> tfvm(new fvMatrix<scalar>(vf,
                                                        vf.dimensions() * dimVol / dimTime));

        const fvMesh &mesh_ = vf.mesh();

        // NO MOVING MESH CORRECTION!
        if(mesh_.moving())
        {
            FatalErrorInFunction<< "Moving mesh not supported for this scheme!" <<exit(FatalError);
        }

        tfvm.ref().diag() = rDeltaTGamma * mesh_.V();
        tfvm.ref().source() = rDeltaTGamma * mesh_.V() * vf.oldTime().primitiveField() + this->getX(vf, currentLoop - 1).ref().primitiveField() * mesh_.V(); // x is from current iteration (currentLoop is +1 to normal loops.)
       
        return tfvm;
    };
    tmp<fvMatrix<vector>>
    SDIRKMethods::ddt(const GeometricField<vector, fvPatchField, volMesh> &vf)
    {

        dimensionedScalar rDeltaTGamma("rDeltaTGamma", dimless / dimTime, 1.0 / (vf.time().deltaT().value() * this->getGamma()));

        initializeRKMethod(vf);
        tmp<fvMatrix<vector>> tfvm(new fvMatrix<vector>(vf,
                                                        vf.dimensions() * dimVol / dimTime));

        const fvMesh &mesh_ = vf.mesh();
        // NO MOVING MESH CORRECTION!
        //check if mesh moved
        if(mesh_.moving())
        {
            FatalErrorInFunction<< "Moving mesh not supported for this scheme!" <<exit(FatalError);
        }

        tfvm.ref().diag() = rDeltaTGamma * mesh_.V();
        tfvm.ref().source() = rDeltaTGamma * mesh_.V() * vf.oldTime().primitiveField() + this->getX(vf, currentLoop - 1).ref().primitiveField() * mesh_.V(); // x is from current iteration (currentLoop is +1 to normal loops.)
        return tfvm;
    };


template <class TYPE>
    void SDIRKMethods::initializeRKMethod(const GeometricField<TYPE, fvPatchField, volMesh> &vf)
    {
        

        bool found = false;
        label vfIndex = -1;

        findFieldPtr(found, vfIndex, vf);

        if (!found)
        {
            // Allocate fList for variable vf.
            allocateF(vfIndex,vf);
        }
        
        //for all loops except the first one.
        if (currentLoop > 1)
        {
            // This would be better served at the end of the RK loop, logically. However, to keep everything in one place, we must do an update based on previous values before solving for C.
            updateFFromPreviousCalc(vfIndex,vf);
        }
    };

    //Find if current Field (such as U) is already stored in this class (as it is not field specific but global definition)
    template<class TYPE>
    void SDIRKMethods::findFieldPtr(bool &found, label &foundIndex, const GeometricField<TYPE, fvPatchField, volMesh> &vf)
    {
         found = false;
         foundIndex = -1;
        forAll(fieldList, i)
        {
            if (fieldList[i] == vf.name()) //non unique names are really bad!
            {
                found = true;
                foundIndex = i;
                break;
            }
        }
    };

     template<class TYPE>
    tmp<GeometricField<TYPE, fvPatchField, volMesh>> SDIRKMethods::getX(const GeometricField<TYPE, fvPatchField, volMesh> &vf,label index)
    {
        
        bool found = false;
        label foundIndex = -1;

        findFieldPtr(found, foundIndex, vf);

        if (!found)
        {
            // Allocate fList for variable vf.
            allocateF(foundIndex, vf);
        }

        if(index > numStages-1)
        {
            FatalErrorInFunction<< "index out of bounds see SDIRKMETHODS.C" <<exit(FatalError);
        }
        
            dimensionedScalar rDeltaT("rDeltaT", dimless / dimTime, 1.0 / (vf.time().deltaT().value()));
            dimensionedScalar myZERO("ZERO",dimless,0.0);

        tmp<GeometricField<TYPE, fvPatchField, volMesh>> tX(new GeometricField<TYPE, fvPatchField, volMesh>(IOobject(
                                                      "X",
                                                      vf.time().timeName(),
                                                      vf.mesh(),
                                                      IOobject::NO_READ,
                                                      IOobject::NO_WRITE),
                                                      vf*myZERO*rDeltaT));

        GeometricField<TYPE, fvPatchField, volMesh> &X = tX.ref();

        updateFFromPreviousCalc(foundIndex, vf);

        
        for (int j = 0; j < index; j++)
        {    
            GeometricField<TYPE, fvPatchField, volMesh>& fj = static_cast<GeometricField<TYPE, fvPatchField, volMesh>&>(fList[foundIndex][j]);
            forAll(X, celli)
            {
                X[celli] += (1.0 / gamma) * this->getAFunction(index, j) * fj[celli];
            }
        }
        
        return tX;
};

template<class TYPE>
    void SDIRKMethods::updateFFromPreviousCalc(label vfIndex ,const GeometricField<TYPE, fvPatchField, volMesh> &vf)
    {

        label prevIndex = currentLoop - 2; // currentLoop Indexes from 1 to inclusive numStages store two times ago


        dimensionedScalar myZERO("0",dimless,0.0);
        if (!fList[vfIndex](prevIndex))
        {
             dimensionedScalar rDeltaTGamma("rDeltaTGamma", dimless / dimTime, 1.0 / (vf.time().deltaT().value() * gamma));
       
           
            fList[vfIndex].set(prevIndex, new  GeometricField<TYPE, fvPatchField, volMesh>(IOobject(
                                                                    vf.name() + "_f" + word(prevIndex),
                                                                    vf.time().timeName(),
                                                                    vf.db(),
                                                                    IOobject::NO_READ,
                                                                    IOobject::NO_WRITE,
                                                                    vf.registerObject()),
                                                                    vf*myZERO*rDeltaTGamma));        
        // allocate values
        //  Find SumR (X)
          
        GeometricField<TYPE, fvPatchField, volMesh>& fj = static_cast<GeometricField<TYPE, fvPatchField, volMesh>&>(fList[vfIndex][prevIndex]);
  
        fj = (vf - vf.oldTime()) * rDeltaTGamma - this->getX(vf, prevIndex).ref(); // -X is from last iteration.
        
        }
    };

 template<class TYPE>
    double SDIRKMethods::findError(GeometricField<TYPE, fvPatchField, volMesh>& vf) //fieldList[fieldi],fList[fieldi]
    {
        bool didFound = false;
        label fieldi = -1;
        findFieldPtr(didFound, fieldi, vf);
        if(!didFound)
        {
            FatalErrorInFunction<< "Did not find field named " << vf.name()  << " see SDIRKMETHODS.C"<< exit(FatalError);
        }
              
            dimensionedScalar rDeltaTGamma("rDeltaTGamma", dimless / dimTime, 1.0 / (vf.time().deltaT().value() * gamma));
            dimensionedScalar myZERO("ZERO",dimless,0.0);
 tmp<GeometricField<TYPE, fvPatchField, volMesh>> tnewf(new GeometricField<TYPE, fvPatchField, volMesh>(IOobject(
                                                      "newf",
                                                      vf.time().timeName(),
                                                      vf.mesh(),
                                                      IOobject::NO_READ,
                                                      IOobject::NO_WRITE),
                                                  vf*rDeltaTGamma*myZERO));
GeometricField<TYPE, fvPatchField, volMesh> & newf = tnewf.ref();
            newf = (vf - vf.oldTime()) * rDeltaTGamma - this->getX(vf, numStages-2).ref(); // -X is from last iteration. numStages -1 should include access to all stored fs.
     tmp<GeometricField<TYPE, fvPatchField, volMesh>> terror(new GeometricField<TYPE, fvPatchField, volMesh>(IOobject(
                                                      "error",
                                                      vf.time().timeName(),
                                                      vf.mesh(),
                                                      IOobject::NO_READ,
                                                      IOobject::NO_WRITE),
                                                  vf*myZERO));        
        GeometricField<TYPE, fvPatchField, volMesh> & error = terror.ref();
    dimensionedScalar bDiff("bdiff",dimless,b[numStages-1]-bHat[numStages-1]);
    error =  (vf.time().deltaT() *bDiff*newf);

            //GeometricField<scalar, fvPatchField, volMesh> errorU = ((vf - vf.oldTime()) * rDeltaTGamma - RKSource) * runTime.deltaT() * timeSchemePtr->getBhatFunction(rkloopCount); // Newest "R"
            for (label bi = 0; bi < numStages-1; bi++)
            {
                GeometricField<TYPE, fvPatchField, volMesh>& fList2 = static_cast<GeometricField<TYPE, fvPatchField, volMesh>&>(fList[fieldi][bi]);
   
                error = error + vf.time().deltaT() *(b[bi]-bHat[bi])*fList2;
            }
            
            //normalize error

            scalarField magField = Foam::mag(error);

            forAll(error,celli)
            { 

              scalar normalizedError = magField[celli]/(mag(vf[celli])+VSMALL);
              scalar cutOff = absoluteRelativeErrorSwitch; //~ switch location
              scalar sharpness = 2;
              scalar alpha = (vf[celli] == TYPE()) ? 0 : 1.0/(Foam::pow(cutOff/(mag(vf[celli])+VSMALL),sharpness)+1.0);
               // setComponent(error[celli],i) = alpha*component(normalizedError,i) + (1.0-alpha)*component(error[celli],i);
                magField[celli] = alpha*normalizedError + (1.0-alpha)*magField[celli];
            }

            double errorScalar = 0;
            if(vf.size() > 0)
            errorScalar = sumSqr(magField)/vf.size(); //eq 105
            //For parallel processing
            reduce(errorScalar, sumOp<scalar>());
            errorScalar = Foam::sqrt(errorScalar); //no sqrt until now.

     

            return errorScalar;
    };

    template<class TYPE>
    void SDIRKMethods::allocateF(label& foundIndex, const GeometricField<TYPE, fvPatchField, volMesh> &vf)
    {
        
        fieldList.append(vf.name());

        // resize arrays to accept new variable.
        fList.resize(fList.size() + 1);                                          // add one to fList
        fList.set(fList.size() - 1, new PtrList<IOobject>(numStages - 1)); // for SDIRK method need stages -1 fs

        foundIndex = fList.size() - 1;
    }
}