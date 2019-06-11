#include "mex.h"
#include "tmwtypes.h"

#if NAN_EQUALS_ZERO
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif

#define O(a,b) orientations[(a-1) + (b-1)*nRows]
#define C(a,b) curvature[(a-1) + (b-1)*nRows]

double fmax(double a, double b){
    return a > b ? a : b;
}
double fmin(double a, double b){
    return a < b ? a : b;
}
double fabs(double a){
    return a > 0.0 ? a : -a;
}

                  
//curvature = CSHRMgetCurvatureMex(orientations)

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
    const mwSize * dimsOrientations;
    uint16_T nRows,nCols;
    uint16_T r,c;
    uint8_T oCount;
    int8_T dr,dc;
    double * curvature,*orientations;
    double oLeft, oRight,dOrLeft,dOrRight;
    
    if (nrhs != 1) {
        mexErrMsgTxt("One input arguments required.");
    } 
    // read input    
    orientations = mxGetPr(prhs[0]);
    dimsOrientations = mxGetDimensions(prhs[0]);
    
    nRows = dimsOrientations[0];
    nCols = dimsOrientations[1];
    
    // create output
    
    plhs[0] = mxCreateDoubleMatrix(nRows,nCols,mxREAL);
   
    
    curvature = mxGetPr(plhs[0]);
    // compute curvature
    
    for(r = 2; r <= nRows-1;r++){        
        for(c = 2;c <= nCols-1;c++){
            if(O(r,c) >= 0.0){
                oCount = 0;
                for(dr = -1;dr < 2;dr++){
                    for(dc = -1;dc < 2;dc++){
                        if(dr == dc && dr == 0) continue;
                        if(O(r+dr,c+dc) > 0.0){
                            if(oCount > 0){
                                oLeft = O(r+dr,c+dc);
                            }else{
                                oRight = O(r+dr,c+dc);
                            }
                            oCount++;
                        }
                    }
                }
                dOrLeft = O(r,c) - oLeft;
                dOrLeft = (dOrLeft > 90.0 ? dOrLeft-180.0 : dOrLeft);
                dOrLeft = (dOrLeft < -90.0 ? dOrLeft+180.0 : dOrLeft);
                
                dOrRight = oRight - O(r,c);
                dOrRight = (dOrRight > 90.0 ? dOrRight-180.0 : dOrRight);
                dOrRight = (dOrRight < -90.0 ? dOrRight+180.0 : dOrRight);
                
                C(r,c) = (oCount == 2 ? fabs(dOrLeft+dOrRight)/2.0 : -1.0);

            }else{
                C(r,c) = -1.0;
            }
        }
    }
}

/*Copyright (c) 2016. Rafael Reisenhofer
  Part of CoShREM Toolbox v1.1
  Built Mon, 11/01/2016
  This is Copyrighted Material*/
