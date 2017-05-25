/* 
[f0] = theory_EIFfpt(params,x0,mu_in,xi,r0,P0,freq)
Calculates the first passage time density of EIF, which can be used to calculate the single-neuron spike train power spectrum.
For derivation of this method, see Richardson Phys. Rev. E 2007.

The SDE is:
C*V' = gL(VL-V)+gL*Delta*exp((V-VT)/Delta)+mu+gL*sigma*sqrt(2*C/gL)*gaussianwhitenoise(t);
with reset at Vr, hard threshold at Vth, lower boundary at Vlb and linear decay of the membrane potential from threshold to reset

freq is a vector of frequencies 
 */

#include "mex.h"
#include "math.h"
#include "complex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

int f, k, Nloop, Nfreq, kre, m;
double sigma2, mu_in, gL, C, Delta, VT, VL, Vth, Vlb, dV, Vr, tref, V, w, pi;
double *params, *freq, *f0r, *f0i;
double df, fmax;

/* get inputs */
params = mxGetPr(prhs[0]);
mu_in = mxGetScalar(prhs[1]);
sigma2 = mxGetScalar(prhs[2]);
freq = mxGetPr(prhs[3]);
Nfreq = mxGetN(prhs[3]);
m = mxGetM(prhs[3]);

if (Nfreq==1 && m!=1) {
    Nfreq = m;
    m = 1;
}

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

Nloop = (int)floor((Vth-Vlb)/dV); /* number of bins */
kre = (int)round((Vr-Vlb)/dV); /* index of reset potential */
pi=3.14159265;

/* allocate outputs */
plhs[0] = mxCreateDoubleMatrix(Nfreq,1,mxCOMPLEX);

f0r = mxGetPr(plhs[0]);
f0i = mxGetPi(plhs[0]);

double *I0, *G, *alpha, *beta;
double _Complex *pf, *po, *jf, *jo;

pf=mxMalloc(Nloop*sizeof(double _Complex));
po=mxMalloc(Nloop*sizeof(double _Complex));
jf=mxMalloc(Nloop*sizeof(double _Complex));
jo=mxMalloc(Nloop*sizeof(double _Complex));
I0=mxMalloc(Nloop*sizeof(double));
G=mxMalloc(Nloop*sizeof(double));
alpha=mxMalloc(Nloop*sizeof(double));
beta=mxMalloc(Nloop*sizeof(double));

/* initialize arrays for integration */

V=Vlb;
for (k=0; k<Nloop; k++){
     
    I0[k] = gL*(VL-V)+mu_in+gL*Delta*exp((V-VT)/Delta);
    G[k] = -I0[k]/(gL*sigma2);
    alpha[k] = exp(dV*G[k]);
    if (G[k]==0) beta[k] = 1/sigma2;
    else beta[k] = (alpha[k]-1)/(G[k]*sigma2);
    
    jf[k] = 0;
    jo[k] = 0;
    pf[k] = 0;
    po[k] = 0;
    
    V += dV;
}


for (f=0; f<Nfreq; f++) {
    
/* initialize flux and probability to 0 */ 
jf[Nloop-1] = 1+0*I; // boundary condition
jo[Nloop-1] = 0+0*I;
pf[Nloop-1] = 0+0*I;
po[Nloop-1] = 0+0*I;
    
w = freq[f]*2*pi;

V = Vth;

for (k = Nloop-2; k>=0; k--){
       
    /* boundary conditions */
    jf[k] = jf[k+1]+dV*I*w*pf[k+1];
    pf[k] = pf[k+1]*alpha[k+1]+C/gL*jf[k+1]*beta[k+1];
    
    /* inhomogenous component - current modulations */
    jo[k] = jo[k+1]+dV*I*w*po[k+1]; if (k==kre) jo[k] -= cexp(-I*w*tref);
    po[k] = po[k+1]*alpha[k+1]+C/gL*jo[k+1]*beta[k+1];
    
    V -= dV;
    
}

f0r[f] = creal(-jo[0]/jf[0]);
f0i[f] = cimag(-jo[0]/jf[0]);

} // frequency loop

mxFree(pf);
mxFree(po);
mxFree(jf);
mxFree(jo);
mxFree(I0);
mxFree(G);
mxFree(alpha);
mxFree(beta);

} //mexFunction