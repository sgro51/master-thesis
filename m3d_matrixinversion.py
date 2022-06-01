#!/usr/bin/env python

import sys
import numpy as np
import scipy
from scipy.sparse import csr_matrix

def read_PETSc_matrix(file):
    f = open(file,'rb')
    id = np.fromfile(f, dtype=">i4", count=1)     # PETSc matrix cookie
    nrow   = int(np.fromfile(f, dtype=">i4", count=1))      # number of rows
    ncol   = int(np.fromfile(f, dtype=">i4", count=1))      # number of columns
    nnzmat, = np.fromfile(f, dtype=">i4", count=1)      # number of nonzeros
    nnzrow = np.fromfile(f, dtype=">i4", count=nrow)   # number of nonzeros per row
    colind = np.fromfile(f, dtype=">i4", count=nnzmat) # column indices
    aij    = np.fromfile(f, dtype=">f8", count=nnzmat) # nonzeros
    f.close()

    rowind=np.zeros(nnzmat)
    start = 0
    for i in range(nrow):
         if nnzrow[i] != 0:
             end = start + nnzrow[i]
             rowind[start:end] = i
             start = end
    A = csr_matrix((aij,(rowind, colind)),shape=(nrow, ncol), dtype=">f8")
    return A

def save_as_PETSc_matrix(output_file_name, A_csr):

    # open file
    f = open(output_file_name, 'wb+')
    # header
    # petsc matrix id 1211216
    np.asarray(1211216, dtype='>i4').tofile(f)
    # number of rows and columns
    nrow, ncol = A_csr.shape
    np.asarray(nrow, dtype='>i4').tofile(f)
    np.asarray(nrow, dtype='>i4').tofile(f)
    # number of nonzeros
    nnzmat = A_csr.getnnz(axis=None)
    np.asarray(nnzmat, dtype='>i4').tofile(f)
    #number of nonzeros per row
    nnzrow=A_csr.getnnz(axis=1)
    np.asarray(nnzrow, dtype='>i4').tofile(f)
    #column indices of nonzeros
    colind = A_csr.indices
    np.asarray(colind, dtype='>i4').tofile(f)
    #nonzeros
    aij = A_csr.data
    np.asarray(aij, dtype='>f8').tofile(f)
    f.close()

def invert_PETSc_matrix(filepathin, filepathout):
    Ai = read_PETSc_matrix(filepathin)
    nnzrow = Ai.getnnz(axis=1)
    Dv_vec = []
    start = 0
    i = 0
    while i < len(nnzrow):
        if nnzrow[i] != 0:
            end = start + nnzrow[i] #size of block
            A_inv = np.linalg.inv(csr_matrix.toarray(Ai[start:end,start:end])) #invert block matrix
            Dv= (-1)*(A_inv - np.identity((end-start))) #compute Dv*dt from Ai_inv=(I-dt*Dv)
            Dv_vec.append(Dv) #creates vector containing all vertical diffusion block matrices
            start = end
        i = i + nnzrow[i]
    Dv_csr=scipy.sparse.block_diag(Dv_vec, format='csr') #converts block matrices into sparse csr matrix again
    save_as_PETSc_matrix(filepathout, Dv_csr)
    return Dv_csr

#
#   main
#
if __name__ == "__main__":
    # usage: m3d_matrixinversion.py filepathin filepathout count
    if len(sys.argv) < 4:
        print("\nNot enough input arguments!")
        print("Usage: m3d_matrixinversion.py filepathin matrixname count \n Example: m3d_matrixinversion.py 1dt/Ai2_%02d.petsc 1dt/Di2_%02d.petsc 12")
        print("Computes the vertical diffusion matrix (multiplied by dt) by inverting implicit transport matrices")
        exit()
    filepathin = sys.argv[1]
    filepathout = sys.argv[2]
    count = int(sys.argv[3])

    for imat in range(0, count):
        filepathin = filepathin % imat
        filepathout = filepathout % imat
        invert_PETSc_matrix(filepathin, filepathout)





