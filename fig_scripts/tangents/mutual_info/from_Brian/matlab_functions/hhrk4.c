/* [tout yout] = hhrk4(stim, tspan, dt, y)

C code that implements a Hodgekin-Huxley model neuron 
via a Runge Kutta fourth order ODE solver.

To be accessed via the following MatLab function:
[tout yout] = hhrk4(stim, tspan, dt)

Inputs:
STIM - stimulus vector that is 2x the length of the output vectors
(and 2x the length of the time vector)
TSPAN - 1x2 vector of inital and final times in msec, e.g. [0 100]
DT - time step in msec
Y - initial values for voltage and gating variables n, m, and h,
	which are commonly [0 0.3177 0.0529 0.5961].

Outputs:
TOUT - 1 x n time vector in ms
YOUT - 4 x n vector of membrane potential (mv) and gating variables n, m, and h

--------------------------------------------------------------------------
Brian Lundstrom 12/04
*/

#include <math.h>
#include <stdio.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double *stim, *tspan, dt, *y, *tout, *yout;
  int nstim, n;
  void hhrk4();
  
  // Check for proper number of arguments.
  if (nrhs != 4) {
    mexErrMsgTxt("Four inputs required.");
  } else if (nlhs != 2) {
    mexErrMsgTxt("Two outputs required.");
  }
   
   // Check to make sure the STIM argument is a 1xn vector.
  if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
	mxGetM(prhs[1]) != 1) {
	mexErrMsgTxt("STIM must be a 1xn vector, e.g. [0 1000].");
}
  
   // Check to make sure the TSPAN argument is a 1x2 vector.
  if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
	mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 2) {
	mexErrMsgTxt("TSPAN must be a 1x2 vector, e.g. [0 1000].");
}
  
  // Check to make sure the DT argument is a scalar.
  if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
	mxGetN(prhs[2])*mxGetM(prhs[2]) != 1) {
	mexErrMsgTxt("DT must be scalar.");
}

// Check to make sure the Y argument is a 1x4 vector.
  if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
	mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != 4) {
	mexErrMsgTxt("Y must be a 1x4 vector, e.g. [0 0.3177 0.0529 0.5961].");
}
  
  // Assign pointers to each input and output.
  stim = mxGetPr(prhs[0]);
  tspan = mxGetPr(prhs[1]);
  dt = mxGetScalar(prhs[2]);
  y = mxGetPr(prhs[3]);
 
  // Check to ensure that STIM length is 2x that of t vector
  double vec_length = ((tspan[1]-tspan[0])/dt+1);
  nstim = mxGetN(prhs[0]);
  if (nstim < (2*vec_length-1)) {
	mexErrMsgTxt("STIM length must be at least 2x of tspan1:dt:tspan2, i.e. STIM sampling is 2/dt.");
}    
  printf("Stimulus time step is %f msec. \n", dt/2.0);
  printf("Output voltage time step is %f msec. \n",dt);
	 	     
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(1, vec_length, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(mxGetN(prhs[3]), vec_length, mxREAL);
  
  tout = mxGetPr(plhs[0]);
  yout = mxGetPr(plhs[1]);
  n = mxGetN(prhs[3]);
  
  /* Call function. */
  hhrk4(stim, tspan, dt, y, tout, yout, n);
}


// hhrk4<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// function HHRK4 gives vectors TOUT and YOUT
//
// NOTE: STIM must be at least 2x length of tspan1:dt:tspan2,
//   i.e. sampling rate of STIM will be 2x that of time
void hhrk4(const double *stim, const double *tspan, const double dt, double *y, double *tout, double *yout, int n)
{

void rungekutta4();
rungekutta4(tspan, y, dt, n, stim, tout, yout);
	
return;
}


// rungekutta4<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// function RUNGEKUTTA4 gives vectors TOUT and YOUT
// the integrator evaluates equations defined in function EQUATIONS
// via a fourth-order Runge-Kutta integration scheme
void rungekutta4(const double *tspan, double *yi, double dt, int n, const double *stim, double *tout, double (*yout)[n])
{
	int i, time_i=1;
	double k1[n], k2[n], k3[n], k4[n], dy[n], tmpy[n], dt2, t;
	
    double y[n];
    for (i=0;i<=(n-1);i++) y[i] = yi[i];
    
	for (i=0;i<=(n-1);i++) yout[0][i] = y[i];		//save initial y vector values
	tout[0] = tspan[0];								//save initial time value
	dt2 = dt/2.0;
	void equations();

	for (t=tspan[0]; t<=tspan[1]; t = t+dt)
	{
		
		/****************************************************************
		implementation of Runge-Kutta Fourth Order routine
		
		dt2 = dt/2.0;  
		k1 = dt * func_ty(t, y)
		k2 = dt * func_ty(t+dt2, y+0.5.*k1)
		k3 = dt * func_ty(t+dt2, y+0.5.*k2)
		k4 = dt * func_ty(t+dt, y+k3)
		ynew = y + 1/6.*(k1 + k4) + 1/3.*(k2 + k3)
		
		Note: MatLab code is used to comment the implementation below.
		****************************************************************/

		equations(t, y, dt, dy, stim);		//~ feval(func_ty, t, y, dt);
		for (i=0;i<=(n-1);i++) {
			k1[i] = dt * dy[i];			//~ k1 = dt*feval(func_ty, t, y, dt);
			}
			
		for (i=0;i<=(n-1);i++) {
			tmpy[i] = y[i]+0.5*k1[i];   //~ y+0.5.*k1
			}
		equations(t+dt2, tmpy, dt, dy, stim); //~ feval(func_ty, t+dt2, y+0.5.*k1, dt);
		for (i=0;i<=(n-1);i++) {
			k2[i] = dt * dy[i];			//~ k2 = dt*feval(func_ty, t+dt2, y+0.5.*k1, dt);
			}
		
		for (i=0;i<=(n-1);i++) {
			tmpy[i] = y[i]+0.5*k2[i];   //~ y+0.5.*k2
			}
		equations(t+dt2, tmpy, dt, dy, stim); //~ feval(func_ty, t+dt2, y+0.5.*k2, dt);
		for (i=0;i<=(n-1);i++) {
			k3[i] = dt * dy[i];			//~ k3 = dt*feval(func_ty, t+dt2, y+0.5.*k2, dt);
			}

		for (i=0;i<=(n-1);i++) {
			tmpy[i] = y[i]+k3[i];		//~ y+k3
			}
		equations(t+dt, tmpy, dt, dy, stim); //~ feval(func_ty, t+dt, y+k3, dt);
		for (i=0;i<=(n-1);i++) {
			k4[i] = dt * dy[i];			//~ k4 = dt*feval(func_ty, t+dt, y+k3, dt);
			}
		
		for (i=0;i<=(n-1);i++) {
			y[i] = y[i] + 1.0/6.0 * (k1[i]+k4[i]) + 1.0/3.0 * (k2[i]+k3[i]);	//substitute new for old y values
			} 

		// create output vectors
		tout[time_i] = t+dt;
		for (i=0;i<=(n-1);i++) yout[time_i][i] = y[i];
		time_i++;
	
	} //end of time_i for loop
return;
}


// equations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// function EQUATIONS outputs value of DY
//
// NOTE: STIM variable is used for time-varying stimuli,
//   where sampling rate of STIM is 2x of 1/dt.
void equations(double t, double *y, double dt, double *dy, double *stim)
{
	double cm, Gl, Gk, Gna, El, Ek, Ena, ts_index;
	// Parameters:
	cm = 10;		// nF / mm2
	Gl = 0.003;  // mS/mm2
	Gk = 0.36;   // mS/mm2
	Gna = 1.2;   // mS/mm2
	El = 10.613-65;// mV
	Ek = -12-65;	// mV 
	Ena = 115-65;	// mV 
   
	//Calculate stimulus value for given time point
	//NOTE: stimulus sampling rate must be 2x that of output sampling rate
	ts_index = rint(2.0/dt*t); //time to stimulus index
	double stimt = stim[(int)ts_index];
		
	//Hodgkin-Huxley differential equations
	dy[0] = -1000*(Gl*(y[0] - El) + Gk*pow(y[1],4)*(y[0] - Ek) + Gna*pow(y[2],3)*y[3]*(y[0] - Ena))/cm + stimt/cm;
	//dy[1] = ((0.1 - 0.01*y[0]) / (exp(1-0.1*y[0]) - 1)) *(1-y[1]) - (0.125 * exp(-.0125*y[0]))*(y[1]);
	//dy[2] = ((2.5 - 0.1*y[0])  / (exp(2.5-0.1*y[0])-1)) *(1-y[2]) - (4 * exp(-y[0]/18))       *(y[2]);
	//dy[3] = (0.07 * exp(-.05*y[0]))                     *(1-y[3]) - (1 / (exp(3-0.1*y[0])+1)) *(y[3]);

    dy[1] = (0.01*(y[0] + 55)) / (1-exp(-.1*(y[0] + 55))) *(1-y[1]) - 0.125 * exp(-.0125*(y[0] + 65))*y[1];
    dy[2] = (0.1*(y[0] + 40)) / (1 - exp(-.1*(y[0] + 40))) *(1-y[2]) - 4 * exp(-.0556*(y[0] + 65))*y[2];
    dy[3] = .07 * exp(-.05*(y[0] + 65)) *(1-y[3]) -  1 / (1 + exp(-.1*(y[0] + 35))) *y[3];
    
return;
}
