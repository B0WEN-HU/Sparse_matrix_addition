/** include file for mex function */
#include "mex.h"
#include "matrix.h"


/** We take sum operation as an example here. You can write a sub or mul operation instead of sum here.*/
/** This function will be called in mexFunction, i.e., the main entry of a matlab mex function. */
void sparse_add(mwIndex* Ap, mwIndex* Bp, mwIndex* Cp, mwIndex* Ai, mwIndex* Bi, mwIndex* Ci, double* Ax, double* Bx, double* Cx, int Cr, int Cc);

/** all mex function should be defined like this*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /** number of rows, columns of matrices A, B , and C */
    int Ar;
    int Ac;
    int Br;
    int Bc;
    int Cr;
    int Cc;
    /** number of non-zeros of A, B and C */
    int Annz, Bnnz, Cnnz;
    
    mwIndex* Ap; /** pointer of matrix A (prhs[0]) */
    mwIndex* Bp; /** pointer of matrix B (prhs[1]) */
    mwIndex* Cp; /** pointer of matrix C (plhs[0], return value) */
    mwIndex* Ai; /** index of matrix A (prhs[0]) */
    mwIndex* Bi; /** index of matrix B (prhs[1]) */
    mwIndex* Ci; /** index of matrix C (plhs[0], return value)*/
    double* Ax; /** non-zeros of matrix A (prhs[0]) */
    double* Bx; /** non-zeros of matrix B (prhs[1]) */
    double* Cx; /** non-zeros of matrix C (plhs[0]) */

  /** check input-output numbers */
    if(nrhs != 2)
        mexErrMsgTxt("Two inputs are required!");
    else if(nlhs != 1)
        mexErrMsgTxt("Only one output is allowed");


    /** check whether the inputs are sparse matrices */
    if(!(mxIsSparse(prhs[0]) && mxIsSparse(prhs[1])))
    {
      mexErrMsgTxt("Two inputs should be sparse matrices");
    }

 
    Ar = mxGetM(prhs[0]); Ac = mxGetN(prhs[0]); Annz = mxGetNzmax(prhs[0]);
    Br = mxGetM(prhs[1]); Bc = mxGetN(prhs[1]); Bnnz = mxGetNzmax(prhs[1]);
    Ap = mxGetJc(prhs[0]); Ai = mxGetIr(prhs[0]); Ax = mxGetPr(prhs[0]);
    Bp = mxGetJc(prhs[1]); Bi = mxGetIr(prhs[1]); Bx = mxGetPr(prhs[1]);

    /** check whether the operation is valid */
    if((Ar != Br) || (Ac != Bc))
      mexErrMsgTxt("The numbers of columns and rows of the two sparse matrices should be matched.");

    /** assign the value of the numbers of rows-cols-non-zeros of C matrix */
    Cr = Ar; Cc = Ac; Cnnz = Annz + Bnnz;
    /** Pointer of C matrix */
    Cp = (mwIndex*)mxMalloc((Cc+1)*sizeof(mwIndex)); 
    Ci = (mwIndex*)mxMalloc((Cnnz)*sizeof(mwIndex));
    Cx = (double*)mxMalloc((Cnnz)*sizeof(double));
    

    /** call the real function. */
    sparse_add(Ap, Bp, Cp, Ai, Bi, Ci, Ax, Bx, Cx, Cr, Cc);
    

    /** clear plhs[0] */
    if(plhs[0] != NULL)
      mxFree(plhs[0]);

    /** create the return sparse matrix from the returned values */
    plhs[0] = mxCreateSparse(Cr, Cc, Cnnz, mxREAL);
    /** set the pointer, index & values of the C matrix */
    mxSetJc(plhs[0], Cp);
    mxSetIr(plhs[0], Ci);
    mxSetPr(plhs[0], Cx);
    

}

/** implementation of the real function. */
/** put your code here. */
/** Cp, Ci, Cx and Cnnz should be computed in this function. */
void sparse_add(mwIndex* Ap, mwIndex* Bp, mwIndex* Cp, mwIndex* Ai, mwIndex* Bi, mwIndex* Ci, double* Ax, double* Bx, double* Cx, int Cr, int Cc)
{
    int i, j;
    /** vectors for scattering */
    /** flags */
    mwIndex* w = (mwIndex*) mxMalloc(Cr*sizeof(mwIndex));
    /** values */
    double* x = (double*) mxMalloc(Cr*sizeof(double));
    int nnz = 0;
    
    /** clear the flags */
  for(i = 0; i < Cr; ++i)
      w[i] = -1;
    
    for(i = 0; i < Cc; ++i)
    {
        Cp[i] = nnz;
        for(j = Ap[i]; j < Ap[i+1]; ++j)
        {
            x[Ai[j]] = Ax[j];
            w[Ai[j]] = i;
            Ci[nnz++] = Ai[j];
        }
        for(j = Bp[i]; j < Bp[i+1]; ++j)
        {
            /** new element */
            if(w[Bi[j]] != i)
            {
                x[Bi[j]] = Bx[j];
                w[Bi[j]] = i;
                Ci[nnz++] = Bi[j];
            }
            else
            {
                x[Bi[j]] += Bx[j];
            }
        }
        
        for(j = Cp[i]; j < nnz; ++j)
        {
            Cx[j] = x[Ci[j]];
        }
    }
    Cp[i] = nnz;
    
    mxFree(w);
    mxFree(x);
}
