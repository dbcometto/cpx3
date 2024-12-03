#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double top, bottom;
    double *inMatrix, *outMatrix;
    size_t ncols, nrows;

    /* check for proper number of arguments */
    if (nrhs != 3)
    {
        mexErrMsgIdAndTxt("MyToolbox:contrastFix:nrhs", "Three inputs required.");
    }
    else if (nlhs != 1)
    {
        mexErrMsgIdAndTxt("MyToolbox:contrastFix:nlhs", "One output required.");
    }

    /* make sure the first input argument is double */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    {
        mexErrMsgIdAndTxt("MyToolbox:contrastFix:notDouble", "Input matrix must be type double.");
    }

    /* make sure the second and third input arguments are doubles */
    if (!mxIsDouble(prhs[1]) ||  mxIsComplex(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1)
    {
        mexErrMsgIdAndTxt("MyToolbox:contrastFix:notScalar", "Top must be a double scalar.");
    }
    else if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1)
    {
        mexErrMsgIdAndTxt("MyToolbox:contrastFix:notScalar", "Top must be a double scalar.");
    }

    /* get the value of the scalar input  */
    top = mxGetScalar(prhs[1]);
    bottom = mxGetScalar(prhs[2]);

/* create a pointer to the real data in the input matrix  */
#if MX_HAS_INTERLEAVED_COMPLEX
    inMatrix = mxGetDoubles(prhs[0]);
#else
    inMatrix = mxGetPr(prhs[0]);
#endif

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[0]);
    nrows = mxGetM(prhs[0]);

    double m = 1.0/(top - bottom);
    double b = 0.0 - m*bottom;

    plhs[0] = plhs[0] = mxCreateDoubleMatrix((mwSize)nrows, (mwSize)ncols, mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    outMatrix = mxGetDoubles(plhs[0]);
#else
    outMatrix = mxGetPr(plhs[0]);
#endif

    for (unsigned short i = 0; i < ncols; i++) {
        for (unsigned short j = 0; j < nrows; j++) {
            unsigned int curr = i*nrows + j;
            if (inMatrix[curr] > top) {
                outMatrix[curr] = 1.0;
            } else if (inMatrix[curr] < bottom) {
                outMatrix[curr] = 0.0;
            } else {
                outMatrix[curr] = m*inMatrix[curr] + b;
            }
        }
    }
}
