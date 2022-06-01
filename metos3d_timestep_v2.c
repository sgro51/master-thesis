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

#include "metos3d_timestep.h"

#undef  kDebugLevel
#define kDebugLevel kDebugLevel2

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepInit"
PetscErrorCode
Metos3DTimeStepInit(Metos3D *metos3d)
{
    PetscFunctionBegin;
    // register event
    PetscLogEventRegister("TimeStepPhi", 0, &metos3d->eventTimeStepPhi);
    // options
    Metos3DUtilOptionsGetScalar(metos3d, "-Metos3DTimeStep", &metos3d->timeStep);
    Metos3DUtilOptionsGetScalar(metos3d, "-Metos3DTimeStepStart", &metos3d->timeStepStart);
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DTimeStepCount", &metos3d->timeStepCount);
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepFinal"
PetscErrorCode
Metos3DTimeStepFinal(Metos3D *metos3d)
{
    PetscFunctionBegin;
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepFinal\n");
    PetscFunctionReturn(0);
}

#undef  kDebugLevel
#define kDebugLevel kDebugLevel3

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepFunction"
PetscErrorCode
Metos3DTimeStepFunction(SNES snes, Vec ynBD, Vec fnBD, void *ctx)
{
    Metos3D     *metos3d    = (Metos3D*)ctx;
    // geometry
    PetscInt    nvec        = metos3d->vectorLength;
    // bgc
    PetscInt    ntracer     = metos3d->tracerCount;
    // load
    PetscInt    nvecloc     = metos3d->vectorLengthLocal;
    // parameter
    PetscInt    nparam      = metos3d->parameterCount;
    PetscReal   *u0         = metos3d->u0;
    // work vars
    PetscInt    itracer;
    Vec         *yin, *yinold, *yout;
    PetscFunctionBegin;
    // create work vectors
    Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &yin, 0.0);
    Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &yinold, 0.0);
    Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &yout, 0.0);
    // fnBD = -Phi(ynBD)+ynBD
    Metos3DUtilVecCopyDiagonalToSeparate(metos3d, ntracer, nvecloc, &ynBD, yin);
    Metos3DUtilVecCopyDiagonalToSeparate(metos3d, ntracer, nvecloc, &ynBD, yinold);
    Metos3DTimeStepPhi(metos3d, yin, yout, nparam, u0);
    for (itracer = 0; itracer < ntracer; itracer++) VecAYPX(yout[itracer], -1.0, yinold[itracer]);
    Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, yout, &fnBD);
    // free work vectors
    VecDestroyVecs(ntracer, &yin);
    VecDestroyVecs(ntracer, &yinold);
    VecDestroyVecs(ntracer, &yout);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepFunction\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepPhi"
PetscErrorCode
Metos3DTimeStepPhi(Metos3D *metos3d, Vec *yin, Vec *yout, PetscInt nparam, PetscReal *u0)
{
    // bgc
    PetscInt    npref       = metos3d->fileFormatPrefixCount;
    PetscInt    ntracer     = metos3d->tracerCount;
    // time step
    PetscScalar t0          = metos3d->timeStepStart;
    PetscScalar dt          = metos3d->timeStep;
    PetscInt    nstep       = metos3d->timeStepCount;
    // work vars
    PetscScalar tj;
    PetscInt    itracer, istep;
    Vec         *ywork;
    PetscFunctionBegin;
    // wait for all processors
    PetscBarrier(PETSC_NULL);
    // start log event
    PetscLogEventBegin(metos3d->eventTimeStepPhi, 0, 0, 0, 0);
    // prepare work vectors
    VecDuplicateVecs(*yin, ntracer, &ywork);
    // initial point in time, project time to [0,1[
    tj = fmod(t0, 1.0);
    // begin bgc, yout not set (yet)
    Metos3DBGCStepBegin(metos3d, tj, dt, yin, yout, nparam, u0);
    // time step loop
    for (istep = 0; istep < nstep; istep++) {
        // compute point in time from deltat and t0
        tj = fmod(t0 + istep*dt, 1.0);
        // work vars
        char filePrefixFormat[PETSC_MAX_PATH_LEN];    
        char filePrefix      [PETSC_MAX_PATH_LEN];
        // file prefix
        if (npref > 1) {
            if ((metos3d->spinupStep + 1)%metos3d->moduloStep[0] == 0) {
                PetscInt modstep = metos3d->moduloStepCount;
                if (modstep > 0) {
                    PetscInt imodstep = metos3d->moduloStep[1];
                    if ((istep+1)%imodstep == 0) {
                        sprintf(filePrefixFormat, "%s%s", metos3d->filePrefix, metos3d->fileFormatPrefix[1]);
                        sprintf(filePrefix, filePrefixFormat, istep);
                        // output
                        Metos3DBGCOutputPrefix(metos3d, filePrefix, ntracer, yin);
                    }
                } else {
                    sprintf(filePrefixFormat, "%s%s", metos3d->filePrefix, metos3d->fileFormatPrefix[1]);
                    sprintf(filePrefix, filePrefixFormat, istep);
                    // output
                    Metos3DBGCOutputPrefix(metos3d, filePrefix, ntracer, yin);
                }
            }
        }
        // yout = PhiStep(yin)
        Metos3DTimeStepPhiStep(metos3d, tj, dt, istep, yin, yout, ywork, nparam, u0);
        for(itracer = 0; itracer < ntracer; itracer++) {
            VecCopy(yout[itracer], yin[itracer]);
        }
    }
    
    // tracer and diag monitor
    if (metos3d->tracerMonitor) Metos3DBGCTracerMonitor(metos3d, yin);
    if (metos3d->diagMonitor) Metos3DBGCDiagMonitor(metos3d);

    // end bgc, yout not set (yet)
    Metos3DBGCStepEnd(metos3d, tj, dt, yin, yout, nparam, u0);
    // free work vector
    VecDestroyVecs(ntracer, &ywork);
    // wait for all processors
    PetscBarrier(PETSC_NULL);
    // stop log event
    PetscLogEventEnd(metos3d->eventTimeStepPhi, 0, 0, 0, 0);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepPhi\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepPhiStep"
PetscErrorCode
Metos3DTimeStepPhiStep(Metos3D *metos3d, PetscScalar t, PetscScalar dt, PetscInt istep, Vec *yin, Vec *yout, Vec *ywork, PetscInt nparam, PetscReal *u0)
{
    PetscBool   flag;
    char        timeint         [PETSC_MAX_PATH_LEN];
    PetscFunctionBegin;
    // time integrator
    sprintf(timeint, "%s", metos3d->timeIntegrator);
    PetscStrcmp("Euler", timeint, &flag);
    if (flag == PETSC_TRUE) {
        // Euler
        Metos3DTimeStepPhiStepEuler(metos3d, t, dt, istep, yin, yout, ywork, nparam, u0);
    }
    else{
        PetscStrcmp("IMEX-2", timeint, &flag);
        if (flag == PETSC_TRUE) {
            // IMEX-2
            Metos3DTimeStepPhiStepIMEX2(metos3d, t, dt, istep, yin, yout, ywork, nparam, u0);
        }
        else{
            PetscStrcmp("IMEX-3", timeint, &flag);
            if (flag == PETSC_TRUE) {
                // IMEX-3
                Metos3DTimeStepPhiStepIMEX3(metos3d, t, dt, istep, yin, yout, ywork, nparam, u0);
            }
            else{
                Metos3DFlag(flag, "This timeintegrator is not implemented (yet)");
            }
        }
    }
    // debug
    Metos3DDebug(metos3d, kDebugLevel, FSSDSE, "Metos3DTimeStepPhiStep", "istep:", istep, "t:", t);
    PetscFunctionReturn(0); 
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepPhiStepEuler"
PetscErrorCode
Metos3DTimeStepPhiStepEuler(Metos3D *metos3d, PetscScalar t, PetscScalar dt, PetscInt istep, Vec *yin, Vec *yout, Vec *ywork, PetscInt nparam, PetscReal *u0)
{   
    // bgc
    PetscInt ntracer        = metos3d->tracerCount;
    // transport
    PetscInt nmat           = metos3d->matrixCount;
    Mat *Ae                 = metos3d->matrixExplicitArray;
    Mat *Ai                 = metos3d->matrixImplicitArray;
    Mat *Aework             = &metos3d->Aework;
    Mat *Aiwork             = &metos3d->Aiwork;
    // work vars
    PetscInt itracer;
    PetscFunctionBegin;
    // bgc
    Metos3DBGCStep(metos3d, t, dt, yin, yout, nparam, u0);
    // transport
    // yout = Ai*(Ae*yin+dt*q(yin))
    Metos3DTransport(metos3d, t, nmat, Ae, ntracer, yin, ywork, Aework);
    for (itracer = 0; itracer < ntracer; itracer++){
        VecAXPY(ywork[itracer], 1.0, yout[itracer]);
    }
    Metos3DTransport(metos3d, t, nmat, Ai, ntracer, ywork, yout, Aiwork);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, FSSDSE, "Metos3DTimeStepPhiStepEuler", "istep:", istep, "t:", t);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepPhiStepIMEX2"
PetscErrorCode
Metos3DTimeStepPhiStepIMEX2(Metos3D *metos3d, PetscScalar t, PetscScalar dt, PetscInt istep, Vec *yin, Vec *yout, Vec *ywork, PetscInt nparam, PetscReal *u0)
{
    //  bgc
    PetscInt ntracer        = metos3d->tracerCount;
    // transport
    PetscInt nmat           = metos3d->matrixCount;
    Mat *Ae2                = metos3d->matrixExplicit2Array;
    Mat *Ai                 = metos3d->matrixImplicitArray;
    Mat *Di2                = metos3d->matrixDiffusionImplicit2Array;
    Mat *Ae2work            = &metos3d->Ae2work;
    Mat *Aiwork             = &metos3d->Aiwork;
    Mat *Di2work            = &metos3d->Di2work;
    // work vars
    PetscInt itracer;
    Vec         *k1_hat_dt;
    PetscFunctionBegin;
    // k1_hat*dt = (Ae2*yin-yin)+dt*q(yin)
    Metos3DBGCStep(metos3d, t, dt, yin, yout, nparam, u0);
    Metos3DTransport(metos3d, t, nmat, Ae2, ntracer, yin, ywork, Ae2work);
    for (itracer = 0; itracer < ntracer; itracer++){
        VecAXPY(ywork[itracer], -1.0, yin[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++){
        VecAXPY(ywork[itracer], 1.0, yout[itracer]);
    }
    VecDuplicateVecs(*ywork, ntracer, &k1_hat_dt); 
    for(itracer = 0; itracer < ntracer; itracer++) {
        VecCopy(ywork[itracer], k1_hat_dt[itracer]);
    }
    // k1*dt = Di2*yin
    Metos3DTransport(metos3d, t-dt, nmat, Di2, ntracer, yin, yout, Di2work);
    // y_hat = ywork = Ai*(yin+dt/2*k1+dt*k1_hat)
    for (itracer = 0; itracer < ntracer; itracer++){
        VecAYPX(yout[itracer], 0.5, yin[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++){
        VecAYPX(yout[itracer],  1.0, k1_hat_dt[itracer]);
    }
    Metos3DTransport(metos3d, t, nmat, Ai, ntracer, yout, ywork, Aiwork);
    // k2_hat*dt = Ae2(t+1)*ywork-ywork+dt*q(ywork)
    Metos3DBGCStep(metos3d, t+dt, dt, ywork, yin, nparam, u0);
    Metos3DTransport(metos3d, t + dt, nmat, Ae2, ntracer, ywork, yout, Ae2work);
    for (itracer = 0; itracer < ntracer; itracer++){
        VecAXPY(yout[itracer], -1.0, ywork[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++){
        VecAXPY(yout[itracer], 1.0, yin[itracer]);
    }
    // yout = yin+dt/2(k1+k2+k1_hat+k2_hat) = y_hat + dt/2*k2_hat-dt/2*k1_hat
    for (itracer = 0; itracer < ntracer; itracer++){
        VecAYPX(yout[itracer], 0.5, ywork[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++){
        VecAXPY(yout[itracer], -0.5, k1_hat_dt[itracer]);
    }
    // free work vector
    VecDestroyVecs(ntracer, &k1_hat_dt);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, FSSDSE, "Metos3DTimeStepPhiStepIMEX2", "istep:", istep, "t:", t);
    PetscFunctionReturn(0);
}


#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepPhiStepIMEX3"
PetscErrorCode
Metos3DTimeStepPhiStepIMEX3(Metos3D *metos3d, PetscScalar t, PetscScalar dt, PetscInt istep, Vec *yin, Vec *yout, Vec *ywork, PetscInt nparam, PetscReal *u0) 
{  
    // bgc
    PetscInt    ntracer     = metos3d->tracerCount;
    // transport
    PetscInt    nmat        = metos3d->matrixCount;
    Mat         *Ae3        = metos3d->matrixExplicit3Array;
    Mat         *Ai3        = metos3d->matrixImplicit3Array;
    Mat         *Ai2        = metos3d->matrixImplicit2Array;
    Mat         *Ae3work    = &metos3d->Ae3work;
    Mat         *Ai3work    = &metos3d->Ai3work;
    Mat         *Ai2work    = &metos3d->Ai2work;
    // work vars
    PetscInt    itracer;
    Vec         *k1_hat_dt;
    Vec         *y_hat;
    Vec         *k2_hat_dt;
    Vec         *y_tilde;
    PetscFunctionBegin;
    // k1_hat_dt = (Ae3*yin-yin)+dt*q(yj)
    Metos3DBGCStep(metos3d, t, dt, yin, yout, nparam, u0);
    Metos3DTransport(metos3d, t, nmat, Ae3, ntracer, yin, ywork, Ae3work);
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(ywork[itracer], -1.0, yin[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(ywork[itracer], 1.0, yout[itracer]);
    }
    VecDuplicateVecs(*ywork, ntracer, &k1_hat_dt);
    for(itracer = 0; itracer < ntracer; itracer++) {
        VecCopy(ywork[itracer], k1_hat_dt[itracer]);
    }
    // y_hat = Ai2(t-1/3dt)*(yin + 2/3dt*k1_hat)
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAYPX(ywork[itracer], 2.0/3.0, yin[itracer]); 
    }
    Metos3DTransport(metos3d, max(t -1.0/3.0*dt, 0.0), nmat, Ai2, ntracer, ywork, yout, Ai2work);
    VecDuplicateVecs(*yout, ntracer, &y_hat);
    for(itracer = 0; itracer < ntracer; itracer++) {
        VecCopy(yout[itracer], y_hat[itracer]);
    }
    // k2_hat*dt = Ae3(t+2/3*dt)*y_hat-y_hat+dt*q(y_hat)
    VecDuplicateVecs(*yout, ntracer, &k2_hat_dt);
    Metos3DBGCStep(metos3d, fmod(t + 2.0/3.0*dt, 1.0), dt, y_hat, k2_hat_dt, nparam, u0);
    Metos3DTransport(metos3d, fmod(t + 2.0/3.0*dt, 1.0), nmat, Ae3, ntracer, y_hat, yout, Ae3work);
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(yout[itracer], -1.0, y_hat[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(k2_hat_dt[itracer], 1.0, yout[itracer]);
    }
    // y_tilde = Ai3*(yin +dt(k2_hat-k1_hat-k2)) = Ai3*(y_hat*(-3/2)+5/2yin+dt*k2_hat)
    for(itracer = 0; itracer < ntracer; itracer++) {
        VecCopy(y_hat[itracer], ywork[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAYPX(ywork[itracer], -3.0/2.0, k2_hat_dt[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(ywork[itracer], +2.5, yin[itracer]);
    }
    Metos3DTransport(metos3d, t, nmat, Ai3, ntracer, ywork, yout, Ai3work);
    VecDuplicateVecs(*yout, ntracer, &y_tilde);
    for(itracer = 0; itracer < ntracer; itracer++) {
        VecCopy(yout[itracer], y_tilde[itracer]);
    }
    // k3_hat*dt = Ae3*y_tilde-y_tilde+dt*q(y_tilde)
    Metos3DBGCStep(metos3d, t, dt, y_tilde, ywork, nparam, u0);
    Metos3DTransport(metos3d, t, nmat, Ae3, ntracer, y_tilde, yout, Ae3work);
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(yout[itracer], -1.0, y_tilde[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(yout[itracer], 1.0, ywork[itracer]);
    }
    // yout = yin + dt*(3/4k2+1/4k3+3/4k2_hat+1/4 k3_hat) = y_tilde*1/4+y_hat*3/2-3/4yin-3/4*dt*k1_hat+1/4*dt*k3_hat+1/2*dt*k2_hat
    VecScale(*yout, 0.25);
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(yout[itracer], 3.0/2.0, y_hat[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(yout[itracer], 1.0/4.0, y_tilde[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(yout[itracer], -3.0/4.0, yin[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(yout[itracer], -3.0/4.0, k1_hat_dt[itracer]);
    }
    for (itracer = 0; itracer < ntracer; itracer++) {
        VecAXPY(yout[itracer], 0.5, k2_hat_dt[itracer]);
    } 
    // free work vectors
    VecDestroyVecs(ntracer, &k1_hat_dt);
    VecDestroyVecs(ntracer, &y_hat);
    VecDestroyVecs(ntracer, &k2_hat_dt);
    VecDestroyVecs(ntracer, &y_tilde);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, FSSDSE, "Metos3DTimeStepPhiStepIMEX3", "istep:", istep, "t:", t);
    PetscFunctionReturn(0);
}
