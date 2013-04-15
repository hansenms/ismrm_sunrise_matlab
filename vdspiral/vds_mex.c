
#include "mex.h" 

#include "vds.c"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

	/*  function [grad] = vds(smax,gmax,Tg,Td,N,F0V,rmax,ngmax) 	*/

 
{
double smax;	/*	Maximum slew rate, G/cm/s		*/
double gmax;	/* 	maximum gradient amplitude, G/cm	*/
double Tgsamp;	/*	Gradient Sample period (s)		*/
double Tdsamp;	/*	Data Sample period (s)			*/
int Nints;	/*	Number of interleaves			*/
double *fieldofview;	/*	FOV coefficients		*/
int Nfov;	/*	Number of FOV coefficients		*/
double krmax;	/*	Maximum k-space extent (/cm)		*/
int ngmax;	/*	Maximum number of gradient samples	*/
double *xgrad;	/* 	X-component of gradient.	*/
double *ygrad;	/*	Y-component of gradient.	*/
double *xgradout;	/* 	X-component of gradient.	*/
double *ygradout;	/*	Y-component of gradient.	*/
double *xgptr, *ygptr;	/* 	temp pointers	*/
int ngrad;	
int count;

/*	Check Arguments 	*/

if (nrhs != 8)
	mexErrMsgTxt("Incorrect Input Argument List");

if (nlhs != 1)
	mexErrMsgTxt("Function has only one output");


/*	Get input parameters	*/

smax = *mxGetPr(prhs[0]);
gmax = *mxGetPr(prhs[1]);
Tgsamp = *mxGetPr(prhs[2]);
Tdsamp = *mxGetPr(prhs[3]);
Nints = (int)(*mxGetPr(prhs[4]));
Nfov = mxGetM(prhs[5])* mxGetN(prhs[5]);
fieldofview = mxGetPr(prhs[5]);
krmax = *mxGetPr(prhs[6]);
ngmax = *mxGetPr(prhs[7]);


/*	Call C-function Here */

calc_vds(smax,gmax,Tgsamp,Tdsamp,Nints,fieldofview,Nfov,krmax,ngmax,
		&xgrad,&ygrad,&ngrad);


/*	Allocate output vector, now that we know size. 	*/

plhs[0] = mxCreateDoubleMatrix(ngrad,1,mxCOMPLEX);
xgradout = mxGetPr(plhs[0]);
ygradout = mxGetPi(plhs[0]);

/*	Copy output into matrix.	*/

bcopy(xgrad,xgradout,ngrad*sizeof(double));
bcopy(ygrad,ygradout,ngrad*sizeof(double));

free(xgrad);
free(ygrad);

}


