/*
 * Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
 * Copyright (C) 2012  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "metos3d_transport.h"

#undef  kDebugLevel
#define kDebugLevel kDebugLevel2

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransportInit"
PetscErrorCode
Metos3DTransportInit(Metos3D *metos3d)
{
    // work vars
    char        transportType[PETSC_MAX_PATH_LEN];    
    PetscBool   flag;
    PetscFunctionBegin;
    // transport type
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DTransportType", transportType);
    PetscStrcmp("Matrix", transportType, &flag);
    if (flag == PETSC_TRUE) {
        // matrix
        Metos3DTransportMatrixInit(metos3d);
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransportInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransportFinal"
PetscErrorCode
Metos3DTransportFinal(Metos3D *metos3d)
{
    // work vars
    char        transportType[PETSC_MAX_PATH_LEN];    
    PetscBool   flag;
    PetscFunctionBegin;
    // transport type
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DTransportType", transportType);
    PetscStrcmp("Matrix", transportType, &flag);
    if (flag == PETSC_TRUE) {
        // matrix
        Metos3DTransportMatrixFinal(metos3d);
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransportFinal\n");
    PetscFunctionReturn(0);
}

#undef  kDebugLevel
#define kDebugLevel kDebugLevel3

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransportMatrixInit"
PetscErrorCode
Metos3DTransportMatrixInit(Metos3D *metos3d)
{
     // work vars
    PetscInt    n_vec     = metos3d->vectorLength;
    PetscInt    n_vec_loc = metos3d->vectorLengthLocal;
    char        inputDirectory  [PETSC_MAX_PATH_LEN];
    char        format          [PETSC_MAX_PATH_LEN];
    char        timeint         [PETSC_MAX_PATH_LEN];
    PetscInt    nmat, imat;
    char        fileFormat      [PETSC_MAX_PATH_LEN];
    PetscBool   flag = PETSC_FALSE; 
    PetscFunctionBegin;
    // input directory
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DMatrixInputDirectory", inputDirectory);
    // matrixCount
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DMatrixCount", &nmat);
    metos3d->matrixCount = nmat;
    // matrix File format
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DMatrixFileFormat", format);
    Metos3DUtilFormatParse(metos3d, format); 
    // time Integrator
    PetscMalloc(PETSC_MAX_PATH_LEN * sizeof(char), &metos3d->timeIntegrator);
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DTimeIntegrator", timeint);
    sprintf(metos3d->timeIntegrator, "%s", timeint);
    PetscStrcmp("Euler", timeint, &flag);
    if (flag == PETSC_TRUE) {
        // matrices for Euler
        Metos3DTransportMatrixLoad(metos3d, "Ae", format, inputDirectory, &metos3d->matrixExplicitArray, &metos3d->Aework);
        Metos3DTransportMatrixLoad(metos3d, "Ai", format, inputDirectory, &metos3d->matrixImplicitArray, &metos3d->Aiwork);
    }
    else{
        PetscStrcmp("IMEX-2", timeint, &flag); 
        if (flag == PETSC_TRUE) {
            // matrices for IMEX-2, Di2=-((Ai2)^(-1)-I)
            Metos3DTransportMatrixLoad(metos3d, "Ae2", format, inputDirectory, &metos3d->matrixExplicit2Array, &metos3d->Ae2work);
            Metos3DTransportMatrixLoad(metos3d, "Ai", format, inputDirectory, &metos3d->matrixImplicitArray, &metos3d->Aiwork);
            Metos3DTransportMatrixLoad(metos3d, "Di2", format, inputDirectory, &metos3d->matrixDiffusionImplicit2Array, &metos3d->Di2work);
        }
        else{ 
            PetscStrcmp("IMEX-3", timeint, &flag); 
            if (flag == PETSC_TRUE) {
                // matrices for IMEX-3
                Metos3DTransportMatrixLoad(metos3d, "Ae3", format, inputDirectory, &metos3d->matrixExplicit3Array, &metos3d->Ae3work);
                Metos3DTransportMatrixLoad(metos3d, "Ai3", format, inputDirectory, &metos3d->matrixImplicit3Array, &metos3d->Ai3work);
                Metos3DTransportMatrixLoad(metos3d, "Ai2", format, inputDirectory, &metos3d->matrixImplicit2Array, &metos3d->Ai2work);
            }
            else{
                Metos3DFlag(flag, "This timeintegrator is not implemented. Please choose between Euler, IMEX-2 and IMEX-3. ");
            } 
        }
    }  
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransportMatrixInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransportMatrixLoad"
PetscErrorCode
Metos3DTransportMatrixLoad(Metos3D *metos3d, char *Aname, char *format, char *inputDirectory, Mat **array, Mat *A_work)
{ 
    // work var 
    PetscInt    n_vec     = metos3d->vectorLength;
    PetscInt    n_vec_loc = metos3d->vectorLengthLocal;
    char        fileFormat  [PETSC_MAX_PATH_LEN];
    char        filePath    [PETSC_MAX_PATH_LEN];
    PetscInt    nmat      = metos3d->matrixCount; 
    PetscInt    imat;
    Mat         *A;
    PetscFunctionBegin;
    // create memory for mat
    PetscMalloc(nmat*sizeof(Mat), &A);
    *array = A; //save in metos3d
    // load all mats
    for (imat = 0; imat < nmat; imat++)
    { // file path
        sprintf(fileFormat, "%s%s%s", "%s", Aname, format); //gives string fileFormat="%sformat"=%sAname_$02d.petsc"
        sprintf(filePath, fileFormat, inputDirectory, imat); //gives string filePath="inputDirectoy/Aname_0imat.petsc"
        // create matrix, set size and type
        MatCreate(metos3d->comm, &A[imat]);
        MatSetSizes(A[imat], n_vec_loc, n_vec_loc, n_vec, n_vec);
        MatSetType(A[imat], MATAIJ);
        // load matrix
        Metos3DUtilMatrixLoad(metos3d, filePath, &A[imat]);
    }
    // work matrix
    MatDuplicate(A[0], MAT_DO_NOT_COPY_VALUES, A_work);
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransportMatrixLoad\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransportMatrixFinal"
PetscErrorCode
Metos3DTransportMatrixFinal(Metos3D *metos3d)
{
    // work vars
    PetscInt    nmat    = metos3d->matrixCount;
    PetscInt    imat;
    PetscBool   flag    = PETSC_FALSE;
    char        timeint   [PETSC_MAX_PATH_LEN]; 
    PetscFunctionBegin;
    // time integrator
    sprintf(timeint, "%s", metos3d->timeIntegrator);
    PetscStrcmp("Euler", timeint, &flag);
    if (flag == PETSC_TRUE) {
        // Euler free matrices
        for (imat = 0; imat < nmat; imat++)
        {   MatDestroy(&metos3d->matrixExplicitArray[imat]);
            MatDestroy(&metos3d->matrixImplicitArray[imat]);
        }
        PetscFree(metos3d->matrixExplicitArray);
        PetscFree(metos3d->matrixImplicitArray);
        // free work matrices
        MatDestroy(&metos3d->Aework);
        MatDestroy(&metos3d->Aiwork);
    }
    else{
        PetscStrcmp("IMEX-2", timeint, &flag);
        if (flag == PETSC_TRUE) {
            // IMEX-2 free matrices
            for (imat = 0; imat < nmat; imat++)
            { 
                MatDestroy(&metos3d->matrixExplicit2Array[imat]);
                MatDestroy(&metos3d->matrixImplicitArray[imat]);
                MatDestroy(&metos3d->matrixDiffusionImplicit2Array[imat]);
            }
            PetscFree(metos3d->matrixExplicit2Array);
            PetscFree(metos3d->matrixImplicitArray);
            PetscFree(metos3d->matrixDiffusionImplicit2Array);
            // free work matrices
            MatDestroy(&metos3d->Ae2work);
            MatDestroy(&metos3d->Aiwork);
            MatDestroy(&metos3d->Di2work);
        }
        else{
            PetscStrcmp("IMEX-3", timeint, &flag);
            if (flag == PETSC_TRUE) {
                // IMEX-3 free matrices
                for (imat = 0; imat < nmat; imat++)
                { 
                    MatDestroy(&metos3d->matrixExplicit3Array[imat]);
                    MatDestroy(&metos3d->matrixImplicit2Array[imat]);
                    MatDestroy(&metos3d->matrixImplicit3Array[imat]);
                }
                PetscFree(metos3d->matrixExplicit3Array);
                PetscFree(metos3d->matrixImplicit2Array);
                PetscFree(metos3d->matrixImplicit3Array);
                // free work matrices
                MatDestroy(&metos3d->Ae3work);
                MatDestroy(&metos3d->Ai2work);
                MatDestroy(&metos3d->Ai3work);
            }
            else{
                Metos3DFlag(flag, "This timeintegrator is not implemented (yet)");
            }
        }
    }
    PetscFree(metos3d->timeIntegrator);
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransportMatrixFinal\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransport"
PetscErrorCode
Metos3DTransport(Metos3D *metos3d, PetscScalar t, PetscInt n_mat, Mat *A, PetscInt n_tracer, Vec *y_in, Vec *y_out, Mat *A_work)
{
    // work vars
    PetscInt    i, i_alpha, i_beta,n,m;
    PetscScalar alpha, beta;
    PetscFunctionBegin;
    // interpolate
    Metos3DUtilInterpolate(metos3d, t, n_mat, &i_alpha, &alpha, &i_beta, &beta);
    // A_work = A[i_alpha]
    // A_work = alpha * A_work
    // A_work = A_work + beta * A[i_beta]
    MatCopy(A[i_alpha], *A_work, SAME_NONZERO_PATTERN);
    MatScale(*A_work, alpha);
    MatAXPY(*A_work, beta, A[i_beta], SAME_NONZERO_PATTERN);
    // apply matrix
    for (i=0; i<n_tracer; i++) {
        MatMult(*A_work, y_in[i], y_out[i]);
    }
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransport\n");
    PetscFunctionReturn(0);
}
