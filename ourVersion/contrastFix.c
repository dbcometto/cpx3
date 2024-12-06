// Install MinGW-w64, then use configuremingw from MATLAB
// Compile with: mex CFLAGS="-Ofast -flto=auto -mfpmath=sse" contrastFix.c

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    register double top, bottom, max;
    static double *inMatrix, *outMatrix;
    size_t ncols, nrows;

    // /* check for proper number of arguments */
    // if (nrhs != 4)
    // {
    //     mexErrMsgIdAndTxt("MyToolbox:contrastFix:nrhs", "Four inputs required.");
    // }
    // else if (nlhs != 1)
    // {
    //     mexErrMsgIdAndTxt("MyToolbox:contrastFix:nlhs", "One output required.");
    // }

    // /* make sure the first input argument is double */
    // if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    // {
    //     mexErrMsgIdAndTxt("MyToolbox:contrastFix:notDouble", "Input matrix must be type double.");
    // }

    // /* make sure the second and third input arguments are doubles */
    // if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1)
    // {
    //     mexErrMsgIdAndTxt("MyToolbox:contrastFix:notScalar", "Top must be a double scalar.");
    // }
    // else if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1)
    // {
    //     mexErrMsgIdAndTxt("MyToolbox:contrastFix:notScalar", "Top must be a double scalar.");
    // }

    /* get the value of the scalar inputs */
    top = mxGetScalar(prhs[1]);
    bottom = mxGetScalar(prhs[2]);
    max = mxGetScalar(prhs[3]);

/* create a pointer to the real data in the input matrix  */
#if MX_HAS_INTERLEAVED_COMPLEX
    inMatrix = mxGetDoubles(prhs[0]);
#else
    inMatrix = mxGetPr(prhs[0]);
#endif

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[0]);
    nrows = mxGetM(prhs[0]);

    top = top * max;
    bottom = bottom * max;
    register double m = max / (top - bottom);
    register double b = -m * bottom;

    plhs[0] = plhs[0] = mxCreateDoubleMatrix((mwSize)nrows, (mwSize)ncols, mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    outMatrix = mxGetDoubles(plhs[0]);
#else
    outMatrix = mxGetPr(plhs[0]);
#endif

    for (unsigned short i = 0; i < ncols * nrows; i++)
    {
        if (inMatrix[i] > top)
        {
            outMatrix[i] = max;
        }
        else if (inMatrix[i] < bottom)
        {
            outMatrix[i] = 0.0;
        }
        else
        {
            outMatrix[i] = m * inMatrix[i] + b;
        }
    }
}
