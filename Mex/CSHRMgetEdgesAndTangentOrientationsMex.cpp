#include "mex.h"
#include "tmwtypes.h"

#if NAN_EQUALS_ZERO
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif

#define CI(a,b,c,d) coeffsImag[(a-1) + (b-1)*nRows + (c-1)*nRows*nCols + (d-1)*nRows*nCols*nOris]
#define CR(a,b,c,d) coeffsReal[(a-1) + (b-1)*nRows + (c-1)*nRows*nCols + (d-1)*nRows*nCols*nOris]
#define PO(a,b) pivotOris[(a-1) + (b-1)*nRows]
#define PS(a,b) pivotScales[(a-1) + (b-1)*nRows]
#define E(a,b) edgeness[(a-1) + (b-1)*nRows]
#define O(a,b) orientations[(a-1) + (b-1)*nRows]

double fmax(double a, double b){
    return a > b ? a : b;
}
double fmin(double a, double b){
    return a < b ? a : b;
}
double fabs(double a){
    return a > 0.0 ? a : -a;
}

                  
//[edges,orientations] = CSHRMgetEdgesAndTangentOrientationsMex(coeffsReal,coeffsImag,pivotOris,pivotScales,shearLevel,minContrast)

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
    const mwSize * dimsCoeffs;
    uint16_T nOris,nScales,nRows,nCols;
    uint16_T r,c,po,s,ps,oConeBorder1,oConeBorder2;
    uint16_T shearLevel;
    uint16_T * pivotOris, * pivotScales;
    double * coeffsReal,* coeffsImag;
    double * angles;
    double * edgeness,*orientations;
    double minContrast;
    double lci,rci,pci,pangle,lpangle,rpangle,rpci,lpci;
    int16_T *linc,*rinc;
    if (nrhs != 6) {
        mexErrMsgTxt("Six input arguments required.");
    } 
    // read input    
    coeffsReal = mxGetPr(prhs[0]);
    coeffsImag = mxGetPr(prhs[1]);
    pivotOris = (uint16_T *)mxGetData(prhs[2]);
    pivotScales = (uint16_T *)mxGetData(prhs[3]);
    shearLevel = mxGetPr(prhs[4])[0];
    minContrast = mxGetPr(prhs[5])[0];

    dimsCoeffs = mxGetDimensions(prhs[0]);
    nRows = dimsCoeffs[0];
    nCols = dimsCoeffs[1];
    nOris = dimsCoeffs[2];
    if(mxGetNumberOfDimensions(prhs[0]) > 3){
        nScales = dimsCoeffs[3];
    }else{
        nScales = 1;
    }
        
    oConeBorder1 = (1<<(shearLevel-2)) + 1;
    oConeBorder2 = nOris/2 - (1<<(shearLevel-2));
    linc = (int16_T *) malloc(nOris*sizeof(int16_T));
    rinc = (int16_T *) malloc(nOris*sizeof(int16_T));
    
    for(r = 0;r<nOris;r++){
        linc[r] = -1;
        rinc[r] = 1;
    }            
    
    linc[0] = nOris - 1; //first element
    linc[oConeBorder1] = -2; //oConeBorder1+1-1
    linc[oConeBorder2] = -2;
    linc[oConeBorder1 + nOris/2] = -2; //oConeBorder1+1-1
    linc[oConeBorder2 + nOris/2] = -2;
    
    rinc[nOris-1] = -nOris+1;
    rinc[oConeBorder1-1] = 2;
    rinc[oConeBorder2-1] = 2;
    rinc[oConeBorder1+nOris/2-1] = 2;
    rinc[oConeBorder2+nOris/2-1] = shearLevel == 2 ? -10 : 2;
    
    
    // create output
    
    plhs[0] = mxCreateDoubleMatrix(nRows,nCols,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nRows,nCols,mxREAL);
    
    
    edgeness = mxGetPr(plhs[0]);
    orientations = mxGetPr(plhs[1]);
    
       
    for(r = 1; r <= nRows;r++){        
        for(c = 1;c <= nCols;c++){
            // compute edgeness

            po = PO(r,c);
            ps = PS(r,c);
            pci = CI(r,c,po,ps);
            
            for(s = 1;s<=nScales;s++){
                E(r,c) = E(r,c) + fmin(CI(r,c,po,s),pci) - CR(r,c,po,s) - minContrast;
            }     
            E(r,c) = fmax(0.0,E(r,c))/(nScales*pci);            
            
            
            
            //compute orientation
            if(E(r,c) > 0.0){
                rpangle = lpangle = pangle = po-1;
                lpci = rpci = pci;
                lci = fmax(CI(r,c,po+linc[po-1],ps),0);
                rci = fmax(CI(r,c,po+rinc[po-1],ps),0);



                if ((po == oConeBorder1) || (po == oConeBorder2) || (po == oConeBorder1+nOris/2) || (po == oConeBorder2+nOris/2)){
                    rpci = fmax(CI(r,c,po+1,ps),0.0);
                    rci = fmin(rci,rpci);
                    rpangle = pangle + 1;
                }
                if ((po == oConeBorder1 + 1) || (po == oConeBorder2 + 1) || (po == oConeBorder1+nOris/2 + 1) || (po == oConeBorder2+nOris/2 + 1)){
                    lpci = fmax(CI(r,c,po-1,ps),0.0);
                    lci = fmin(lci,lpci);
                    lpangle = pangle - 1;
                }
                if (rci/rpci > lci/lpci){
                    lci = rpci*lci/lpci;
                    pangle = rpangle;
                    pci = rpci;
                }else{
                    rci = lpci*rci/rpci;
                    pangle = lpangle;
                    pci = lpci;
                }
            
    
                O(r,c) = pangle + (rci - lci)/(4*pci - 2*(rci + lci));
                //O(r,c) = pangle + (rci - lci)/(2*(pci - fmin(rci,lci)));
                if(O(r,c) < 0.0) O(r,c) = O(r,c) + nOris;
                if(O(r,c) >= nOris) O(r,c) = O(r,c) - nOris;
                if(O(r,c) >= (nOris)/2.0) O(r,c) = O(r,c) - nOris/2.0;
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
