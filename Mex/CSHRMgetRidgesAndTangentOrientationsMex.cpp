#include "mex.h"
#include "tmwtypes.h"

#if NAN_EQUALS_ZERO
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif

#define CR(a,b,c,d) coeffsReal[(a-1) + (b-1)*nRows + (2*(c-1)+hasEvenLength)*nRows*nCols + (d-1)*nRows*nCols*2*nOris]
#define CI(a,b,c,d) coeffsImag[(a-1) + (b-1)*nRows + (2*(c-1)+hasEvenLength)*nRows*nCols + (d-1)*nRows*nCols*2*nOris]
#define PW(a) positiveWidths[2*(a-1)+hasEvenLength]
#define EC(a,b) expectedCoeffs[(a-1)+(b-1)*maxWidth] //usage: EC(width,scale)
#define PO(a,b) pivotOris[(a-1) + (b-1)*nRows]
#define PS(a,b) pivotScales[(a-1) + (b-1)*nRows]
#define L(a,b) ridgeness[(a-1) + (b-1)*nRows]
#define O(a,b) orientations[(a-1) + (b-1)*nRows]
#define W(a,b) widths[(a-1) + (b-1)*nRows]
#define H(a,b) heights[(a-1) + (b-1)*nRows]


double fmax(double a, double b){
    return a > b ? a : b;
}
double fmin(double a, double b){
    return a < b ? a : b;
}
double fabs(double a){
    return a > 0.0 ? a : -a;
}
//[edges,orientations] = CSHRMgetRidgesMex(coeffsReal,coeffsImag,pivotOris,pivotScales,positiveWidths,expectedCoeffs,shearLevel,minContrast)

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
    const mwSize * dimsCoeffs, * dimsExpectedCoeffs;
    uint8_T hasEvenLength = 0;
    uint16_T nOris,nScales,nRows,nCols,maxWidth;
    uint16_T r,c,po,s,ps,oConeBorder1,oConeBorder2;
    uint16_T shearLevel;
    uint16_T * pivotOris,*pivotScales,*positiveWidths;
    double * coeffsReal,* coeffsImag,* expectedCoeffs;
    double *ridgeness,*orientations,*widths,*heights;
    double minContrast;
    double lcr,rcr,pcr,pangle,lpangle,rpangle,rpcr,lpcr,normalizationFactor;
    int16_T *linc,*rinc;
    if (nrhs != 8) {
        mexErrMsgTxt("Eight input arguments required.");
    } 
    // read input    
    
    coeffsReal = mxGetPr(prhs[0]);
    coeffsImag = mxGetPr(prhs[1]);
    pivotOris = (uint16_T *)mxGetData(prhs[2]);    
    pivotScales = (uint16_T *)mxGetData(prhs[3]);
    positiveWidths = (uint16_T *)mxGetData(prhs[4]);
    expectedCoeffs = mxGetPr(prhs[5]);
    shearLevel = (uint16_T) mxGetPr(prhs[6])[0];
    minContrast = mxGetPr(prhs[7])[0];

    dimsCoeffs = mxGetDimensions(prhs[0]);
    nRows = dimsCoeffs[0];
    nCols = dimsCoeffs[1];
    nOris = dimsCoeffs[2]/2;
    if(mxGetNumberOfDimensions(prhs[0]) > 3){
        nScales = dimsCoeffs[3];
    }else{
        nScales = 1;
    }
    
    if(nScales == 0) nScales = 1;
    
    dimsExpectedCoeffs = mxGetDimensions(prhs[5]);
    maxWidth = dimsExpectedCoeffs[0];
        
    oConeBorder1 = (1<<(shearLevel-2)) + 1;
    oConeBorder2 = nOris - (1<<(shearLevel-2));
    linc = (int16_T *) malloc(nOris*sizeof(int16_T));
    rinc = (int16_T *) malloc(nOris*sizeof(int16_T));
    
    for(r = 0;r<nOris;r++){
        linc[r] = -1;
        rinc[r] = 1;
    }            
    
    linc[0] = nOris - 1; //first element
    linc[oConeBorder1] = -2; //oConeBorder1+1-1
    linc[oConeBorder2] = -2;
    
    rinc[nOris-1] = -nOris+1;
    rinc[oConeBorder1-1] = 2;
    rinc[oConeBorder2-1] = shearLevel == 2 ? -4 : 2;
    
   
    // create output
    
    plhs[0] = mxCreateDoubleMatrix(nRows,nCols,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nRows,nCols,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nRows,nCols,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nRows,nCols,mxREAL);
    
    ridgeness = mxGetPr(plhs[0]);
    orientations = mxGetPr(plhs[1]);
    widths = mxGetPr(plhs[2]);
    heights = mxGetPr(plhs[3]);
    
    // compute ridgeness    

    for(r = 1; r <= nRows;r++){        
        for(c = 1;c <= nCols;c++){
            po = PO(r,c);
            hasEvenLength = (po-1) % 2;
            po = (po-1)/2+1;
            
            ps = PS(r,c);
            
            W(r,c) = PW(ps);
            H(r,c) = CR(r,c,po,ps)/EC(PW(ps),ps);
            
            normalizationFactor = 0.0;
            for(s = 1;s<=nScales;s++){
                L(r,c) += CR(r,c,po,s);
                normalizationFactor += fmax(fabs(CR(r,c,po,s)),fabs(H(r,c)*EC(PW(ps),s)));
            }              
            L(r,c) = fabs(L(r,c));
            for(s = 1;s<=nScales;s++){
                L(r,c) = L(r,c) - CI(r,c,po,s) - minContrast*fabs(EC(PW(ps),s));
            }
            L(r,c) = fmax(0.0,L(r,c))/normalizationFactor;   

            

            //compute orientation
            if(L(r,c) > 0.0){
                pcr = fabs(CR(r,c,po,ps));  
                lpcr = rpcr = pcr;

                lcr = fabs(CR(r,c,po+linc[po-1],ps));
                rcr = fabs(CR(r,c,po+rinc[po-1],ps));
                
                rpangle = lpangle = pangle = po-1;



               if ((po == oConeBorder1) || (po == oConeBorder2)){
                    rpcr = fmax(CR(r,c,po+1,ps),0.0);
                    rcr = fmin(rcr,rpcr);
                    rpangle = pangle + 1;
                }
                if ((po == oConeBorder1+1) || (po == oConeBorder2+1)){
                    lpcr = fmax(CR(r,c,po-1,ps),0.0);
                    lcr = fmin(lcr,lpcr);
                    lpangle = pangle - 1;
                }
                if (rcr/rpcr > lcr/lpcr){
                    lcr = rpcr*lcr/lpcr;
                    pangle = rpangle;
                    pcr = rpcr;
                }else{
                    rcr = lpcr*rcr/rpcr;
                    pangle = lpangle;
                    pcr = lpcr;
                }
            
    
                //O(r,c) = pangle + (rcr - lcr)/(4*pcr - 2*(rcr + lcr));
                O(r,c) = pangle + (rcr - lcr)/(2*(pcr - fmin(rcr,lcr)));
                if(O(r,c) < 0.0) O(r,c) = O(r,c) + nOris;
                if(O(r,c) >= nOris) O(r,c) = O(r,c) - nOris;
                O(r,c) = O(r,c) + 1;
            } else {
                O(r,c) = -1;
            }            
        }
    }
    delete [] rinc;
    delete [] linc;
}


/*Copyright (c) 2016. Rafael Reisenhofer
  Part of CoShREM Toolbox v1.1
  Built Mon, 11/01/2016
  This is Copyrighted Material*/
