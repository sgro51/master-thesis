# master-thesis

These programs can be seen as an addition to Metos3D (https://github.com/metos3d), the toolkit for marine ecosystem modelling. They add the possibility to choose between different time integrators, namely the semi-implicit Euler scheme, a semi-implicit Runge-Kutta scheme of order 2 (IMEX-2) and a semi-implicit Runge-Kutta-scheme of order 3 (IMEX-3). 

The integration into Metos3D works as follows:
After downloading Metos3D as on https://metos3d.github.io/, one can replace the programs metos3d_type.h, metos3d_transport.h, metos3d_transport.c, metos3d_timestep.h, metos3d_timestep.c (in the file metos3d/simpack/src) by the respective programs metos3d_type_v2.h, metos3d_transport_v2.h, metos3d_transport_v2.c, metos3d_timestep_v2.h, metos3d_timestep_2.c. 

When running the programm, an option file has to be provided. Examples for this can be found in metos3d/model/model/N/option/test.N.option.txt. The option file has to be changed to represent the time integrator. A new example using the IMEX-3 scheme for the N model can therefore be found in option_IMEX-3.txt.

To use the IMEX-2 or IMEX-3 scheme, inverse transport matrices (which correspond to the vertical diffusion matrix) and matrices in different step sizes have to be provided. The matrices with bigger step sizes can be computed with the matrix command. Its usage is

metos3d [-v] matrix [exp|imp] [count] [factor] [file-format-in] [file-format-out].

To also compute the inverse matrices, replace the program metos3d_mod.py in metos3d/metos3d by metos3d_mod_v2.py. Then, inverse matrices can be computed with the matrixinversion command. Its usage is 

metos3d [-v] matrixinversion [count] [file-format-in] [file-format-out].

For the IMEX-2 scheme with a step size Delta t, the matrices for Delta t, the implicit matrix for Delta t /2 and the inverse implicit matrix for Delta t have to be provided.
An example would be the step size 2dt (smallest possible step size for IMEX-2): 
The input directory is specified in the option file by 

-Metos3DMatrixInputDirectory data/TMM/2.8/Transport/Matrix5_4/1dt/.

In the directory metos3d/data/TMM/2.8/Transport/Matrix5_4, create the new matrices by:

metos3d matrix exp 12 2 1dt/Ae_%02d.petsc 1dt/Ae2_%02d.petsc\\
metos3d matrix imp 12 2 1dt/Ai_%02d.petsc 1dt/Ai2_%02d.petsc

metos3d matrixinversion 12 1dt/Ai2_%02d.petsc 1dt/Di2_%02d.petsc.

For the IMEX-3 scheme, the smallest possible step size is 3dt. In this case, implicit and explicit matrices for 3dt and the implicit matrix for 2dt is needed. This is done by:

metos3d matrix exp 12 3 1dt/Ae_%02d.petsc 1dt/Ae3_%02d.petsc

metos3d matrix imp 12 3 1dt/Ai_%02d.petsc 1dt/Ai3_%02d.petsc

metos3d matrix imp 12 2 1dt/Ai_%02d.petsc 1dt/Ai2_%02d.petsc.

The python script m3d_matrixinversion offers a second possibility of creating the vertical diffusion matrices.

The step sizes provided in the option file have to be adapted to the scheme. The smallest possible step size for IMEX-2 is 2dt, which corresponds to the options

 -Metos3DTimeStepCount                              1440
 
 -Metos3DTimeStep                                   0.000694444444444
 
 The smallest step size for IMEX-3 is 3dt, which is displayed in the option file by 
 
 -Metos3DTimeStepCount                              960
 
 -Metos3DTimeStep                                   0.00104166666666.
 
 Bigger step sizes can of course be chosen, but the required matrices have to be prepared accordingly.
