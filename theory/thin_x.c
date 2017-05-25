/* 
[P0,p0,J0,r0,x0] = thin_x(params,x0_out(2),mu_in,xi)
Calculates the stationary density and firing rate of an EIF with a voltage-activated conductance and white noise input
The SDE is:
C*V' = gL(VL-V)+gx*x*(Vx-V)+gL*Delta*exp((V-VT)/Delta)+mu+gL*sigma*sqrt(2*C/gL)*xi(t);
with reset at Vr, hard threshold at Vth, lower boundary at Vlb and linear decay of the membrane potential from threshold to reset
*/

#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

int k,Nloop,kre;
double x0_in, mu_in;
double sigma2, gL, C, Delta, VT, VL, Vth, Vlb, dV, Vr, tref, tau_x, Vx, gx, V;
double *params, *xi, *P0, *p0, *J0, *r0, *x0, *Psp,*j0,*I0,*G,*alpha,*beta;
double p0sum, P0sum;
        

/* get inputs */
params = mxGetPr(prhs[0]);
x0_in = mxGetScalar(prhs[1]);
mu_in = mxGetScalar(prhs[2]);
sigma2 = mxGetScalar(prhs[3]);
xi = mxGetPr(prhs[4]);

gL = params[0];
C = params[1];
Delta = params[2];
VT = params[3];
VL = params[4];
Vth = params[5];
Vlb = params[6];
dV = params[7];
Vr = params[8];
tref = params[9];
tau_x = params[10];
Vx = params[11];
gx = params[12];

Nloop = (int)floor((Vth-Vlb)/dV);  /* number of bins */
kre = (int)round((Vr-Vlb)/dV); /* index of reset potential */


Psp=mxMalloc(Nloop*sizeof(double));
j0=mxMalloc(Nloop*sizeof(double));
I0=mxMalloc(Nloop*sizeof(double));
G=mxMalloc(Nloop*sizeof(double));
alpha=mxMalloc(Nloop*sizeof(double));
beta=mxMalloc(Nloop*sizeof(double));

Psp[Nloop-1] = 0;
j0[Nloop-1] = 1;

/*
double Psp[Nloop]; Psp[Nloop-1] = 0;
double j0[Nloop]; j0[Nloop-1] = 1; // boundary condition
double I0[Nloop];
double G[Nloop];
double alpha[Nloop];
double beta[Nloop];
*/

/* allocate output */
plhs[0] = mxCreateDoubleMatrix(Nloop,1,mxREAL);
plhs[1] = mxCreateDoubleMatrix(Nloop,1,mxREAL);
plhs[2] = mxCreateDoubleMatrix(Nloop,1,mxREAL);
plhs[3] = mxCreateDoubleScalar(0);
plhs[4] = mxCreateDoubleScalar(0);
        
P0 = mxGetPr(plhs[0]);
p0 = mxGetPr(plhs[1]); p0[Nloop-1] = 0;
J0 = mxGetPr(plhs[2]);
r0 = mxGetPr(plhs[3]);
x0 = mxGetPr(plhs[4]);

/* integrate backwards from threshold */
V = Vth;
I0[Nloop-1] = gL*(VL-V)+gx*x0_in*(Vx-V)+mu_in+gL*Delta*exp((V-VT)/Delta);    
G[Nloop-1] = -I0[Nloop-1]/(gL*sigma2);
alpha[Nloop-1] = exp(dV*G[Nloop-1]);

// if (G[Nloop-1]<.000000000001) beta[Nloop-1] = dV/sigma2;
if (G[Nloop-1]==0) beta[Nloop-1] = dV/sigma2;
else beta[Nloop-1] = (alpha[Nloop-1]-1)/(G[Nloop-1]*sigma2);


p0sum = 0;

for (k = Nloop-2; k>0; k--){
    
    V = V-dV;
    
    I0[k] = gL*(VL-V)+gx*x0_in*(Vx-V)+mu_in+gL*Delta*exp((V-VT)/Delta);
    G[k] = -I0[k]/(gL*sigma2);
    alpha[k] = exp(dV*G[k]);
//     if (G[k]<.000000000001) {beta[k] = dV/sigma2; mexPrintf("\naa=%f %f\n", dV/sigma2, (alpha[k]-1)/(G[k]*sigma2)); }
    if (G[k]==0) beta[k] = dV/sigma2;
    else beta[k] = (alpha[k]-1)/(G[k]*sigma2);
    
    
    j0[k] = j0[k+1];
    if (k==kre) j0[k] += -1;
    
    p0[k] = p0[k+1]*alpha[k+1]+C/gL*j0[k+1]*beta[k+1];
    
    p0sum += p0[k];
    Psp[k] = 0;
    
}

/* normalize and compute firing rate and mean gx activation */

(*r0)=1/(dV*p0sum+tref);

P0sum = 0;
(*x0) = 0;
for (k=0; k<Nloop; k++) {
    J0[k] = (*r0)*j0[k];
    if (k >= kre) Psp[k] = (*r0)*tref/(Vth-Vr);
    
    P0[k] = (p0[k]*(*r0))+Psp[k]; // linear spike shape
    P0sum += P0[k];
    (*x0) += P0[k]*xi[k]/tau_x;
}

(*x0) = dV*(*x0)/(dV*P0sum/tau_x);

mxFree(Psp);
mxFree(j0);
mxFree(I0);
mxFree(G);
mxFree(alpha);
mxFree(beta);


}