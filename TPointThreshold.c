/* TPointThreshold.c
*   threshold of an unimodal histogram
*
*         Created:  may/2008   Nicolas COUDRAY and Jean-Luc BUESSLER
*         Last modification:  -
*         Version 0.5
*
*         Copyright: 2008 %AnIM toolbox for HT3DEM   www.ht3dem.org
*         University of Haute Alsace - MIPS Laboratory - www.trop.mips.uha.fr
*
*     reference:
*     N. Coudray, J.-L. Buessler, J.-P. Urban
*     "Robust threshold estimation for images with unimodal histograms"
*     to be publised in "Pattern Recognition Letters" (2010)
*/
#include "mex.h"

void T_Point (double ptX[],double ptY[],double* ptTh,double ptEc[],int nPoints )
{
    int i,indMax=0;
    double coefA, cov, xVar, er;
    double valMax=1e300;
 
    /* variables evaluated recursively */ 
    int nval=0;
    double sumX=0;
    double sumY=0;
    double sumXY=0;
    double sumXX=0;
    double sumYY=0;
    
    /* initialisation of first values */
    ptEc[0]=0;
    ptEc[1]=0;
    
    /* Evaluation of Ec1 error function for points from 1 to k */
    for( i=0;i<nPoints;i++)
    {
        nval++;                     // number of values
        sumX= ptX[i]+sumX;          // sum of X 
        sumY= ptY[i]+sumY;          // sum of Y
        sumXY=ptX[i]*ptY[i]+sumXY;  // sum of products X et Y
        sumXX= ptX[i]*ptX[i]+sumXX;     // sum of squared X
        sumYY=ptY[i]*ptY[i]+sumYY;      // sum of squared Y
        
        if (nval>2)
        {
            xVar=sumXX-(sumX*sumX/nval);// variance de X
            if (xVar ==0)
                mexErrMsgTxt("Error vector X: identical values");
            cov=sumXY-(sumX*sumY/nval); //covariance
            coefA=cov/xVar;             //slope estimation
            ptEc[i]=sumYY - (sumY*sumY )/nval-(coefA * cov);  //Ec1 : sum of squared residuals
       }
    }
    
    /* reinitialization of recursive variables */
    nval=0;
    sumX=0;
    sumY=0;
    sumXY=0;
    sumXX=0;
    sumYY=0;
        
    /* Evaluation of Ec2 error function for points from k+1 to n */
    for(i=nPoints-1;i>0;i--)
    {
        nval++;
        sumX= ptX[i]+sumX;          // sum of X 
        sumY= ptY[i]+sumY;          // sum of Y
        sumXY=ptX[i]*ptY[i]+sumXY;  // sum of products X et Y
        sumXX= ptX[i]*ptX[i]+sumXX;     // sum of squared X
        sumYY=ptY[i]*ptY[i]+sumYY;      // sum of squared Y
        
        if (nval>2)
        {
            xVar=sumXX - (sumX*sumX/nval);// variance de X
            if (xVar ==0)
                mexErrMsgTxt("Error vector X: identical values");
            cov=sumXY-(sumX*sumY/nval); // covariance
            coefA=cov/xVar;             // slope estimation
            er=sumYY -(sumY*sumY )/nval-(coefA * cov);  //sum of squared residuals
            
            ptEc[i-1]= ptEc[i-1] + er ;
        }
        if (valMax>ptEc[i])
        {
            valMax=ptEc[i];
            indMax=i;
        }
   }
    *ptTh = indMax+1;
}

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    double *ptX,*ptY,*ptTh,*ptEc;
    int nX,nY,i;
    int mrowsX,ncolsX,mrowsY,ncolsY;
    mxArray * EcArray;
    
    /* Check input arguments */
    if(nrhs== 0 || nrhs>2) 
        mexErrMsgTxt("One or Two inputs required.");
    
  /*  get the dimensions of the matrix input x and y */
    mrowsY = mxGetM(prhs[0]);
    ncolsY = mxGetN(prhs[0]);
    nY=ncolsY*mrowsY;
    ptY = mxGetPr(prhs[0]);
    
     /* Check sizes */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || ((mrowsY > 1 )&&( ncolsY > 1)))
        mexErrMsgTxt("Input must be vector of doubles.");
        
    if (nrhs==1)
    {
        ptX=mxCalloc(nY,sizeof(double));
        for ( i=0;i<nY;i++)
            ptX[i]=i;
    }
    else
    {
        mrowsX = mxGetM(prhs[1]);
        ncolsX = mxGetN(prhs[1]);
        nX=ncolsX*mrowsX;
        ptX = mxGetPr(prhs[1]);
        
        /* Check sizes */
        if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || ((mrowsX > 1 )&&( ncolsX > 1)))
            mexErrMsgTxt("Input must be vectors of doubles.");
        if(nY!=nX)
            mexErrMsgTxt("Parameter vectors must be of same dimension.");
    }


  /* Create matrix for the return argument. */
    plhs[0] = mxCreateDoubleMatrix(1, 1,mxREAL);
    
    EcArray= mxCreateDoubleMatrix(nX, 1,mxREAL);
    if (nlhs >1)      
        plhs[1]=EcArray;

  /*  create a C pointer to a copy of the output matrix */
    /*pointeur sur les sorties*/
    ptTh = mxGetPr(plhs[0]);
    ptEc = mxGetPr(EcArray);
       /*  appels de la routine en c */
    T_Point(ptX,ptY,ptTh,ptEc,nY);
      
}


