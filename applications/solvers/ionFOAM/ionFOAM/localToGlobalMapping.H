 //Internal faces
            forAll(phi,fluidFacei)
            {
                 label globalFacei = faceAddressList[fluidFacei]-1; //Minus 1!!!
                 phiGlobal[globalFacei] = phi[fluidFacei];
            }

            forAll(meshFluid.boundaryMesh(),fluidBoundaryi)
            {
                fvsPatchScalarField& bci = phi.boundaryFieldRef()[fluidBoundaryi]; 
                label startIndex = bci.patch().start();
                label sizeOfPatch = bci.patch().size(); 

                //If not -1, then boundary exists in global mesh. 
                //Otherwise, fluid boundary data needs to be mapped to internal faces of phiGlobal.
                if(boundaryRegionAddressList[fluidBoundaryi] > -1)
                {
                       fvsPatchScalarField& bcGlobali = phiGlobal.boundaryFieldRef()[boundaryRegionAddressList[fluidBoundaryi]];
                       
                       for(label fluidFacei = startIndex; fluidFacei < startIndex+sizeOfPatch; fluidFacei++)
                        {  
                            label localPatchFacei = fluidFacei-startIndex;
                            label globalFacei = faceAddressList[fluidFacei]; // need to determine sign before flip!

                            globalFacei = globalFacei >= 0 ? globalFacei - 1: -1*globalFacei -1;
                            label globalPatchFacei = globalFacei - bcGlobali.patch().start();
                            bcGlobali[globalPatchFacei] = bci[localPatchFacei];
                        } 
                        
                }
                else
                {
                       for(label fluidFacei = startIndex; fluidFacei < startIndex+sizeOfPatch; fluidFacei++)
                        {  
                            
                            label localPatchFacei = fluidFacei-startIndex;
                            label globalFacei = faceAddressList[fluidFacei]; // need to determine sign before flip!

                            globalFacei = globalFacei >= 0 ? globalFacei - 1: -1*globalFacei -1;
                            phiGlobal[globalFacei] = bci[localPatchFacei];
                        }
                }
            }