/*==========================================================
 *
 * weightCalculator.c calculates sequence weights for MSA
 *
 * Input: MSA in numerical form, matrix dimensions, threshold
 *
 * Output: vector of sequence weights
 *
 * This is a MEX-file for MATLAB. Compile with "mex weightCalculator.c"
 *
 *========================================================*/

#include "mex.h"

/* The computational routine */
void weightCalculator(double *align, double theta, double m, double length, double *W)
{
    mwSize dist;
    int i,j,k,a,b;
    
    /* initialize W to one */
    
    for (i=0; i<m; i++)
        W[i] = 1.0;
    
    /* count seqs below theta dist */
    
    for (i=0; i<m-1; i++) 
        for (j=i+1; j<m; j++)
            {   
                dist = 0.0;
                for (k=0;k<length;k++)
                {
                    a = k+i*length;
                    b = k+j*length;
                    dist += ( align[a] != align[b] );
                }
                if ( 1.0*dist < theta * length )
                {
                    W[i] += 1.0;
                    W[j] += 1.0;
                }
               
            }
    
    /* weight = 1/number of sequences */
   
    for (i=0; i<m; i++)
        W[i] = 1.0/W[i];
    
    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *align;                  /* alignment */
    double theta;                   /* scalar threshold */
    mwSize m;                       /* number of sequences */
    mwSize length;                  /* number of AAs = length of seq */
    double *outMatrix;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Four inputs required: align, theta, M, N.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required: W.");
    }
    
    /* get the value of the alignment  */
    align = mxGetPr(prhs[0]);

    /* get the threshold value theta */
    theta = mxGetScalar(prhs[1]);
    
    /* get dimensions of the input matrix */
    m = mxGetScalar(prhs[2]);
    length = mxGetScalar(prhs[3]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,m,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    weightCalculator(align,theta, m,length,outMatrix);
}
