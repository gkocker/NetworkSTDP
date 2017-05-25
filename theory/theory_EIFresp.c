/* 
[A,x1,V1] = theory_EIFresp(params,x0,mu_in,xi,r0,P0,freq)
Calculates the linear response an EIF with a voltage-activated conductance to a small periodic perturbation in the mean input current
For derivation of this method, see Richardson Phys. Rev. E 2009, which extends the threshold integration method of Richardson Phys. Rev. E 2007 to models
with slow voltage- or calcium-activated conductances.

Since we have set the strength of that additional conductance to 0, this is equivalent to the method of Richardson Phys. Rev. E 2007.

The SDE is:
C*V' = gL(VL-V)+gx*x*(Vx-V)+gL*Delta*exp((V-VT)/Delta)+mu+gL*sigma*sqrt(2*C/gL)*gaussianwhitenoise(t);
with reset at Vr, hard threshold at Vth, lower boundary at Vlb and linear decay of the membrane potential from threshold to reset

freq is a vector of frequencies 
A is the linear response function, x1 is the first-order response of the activation variable for the voltage-activated conductance and V1 is the linear response of the membrane potential

 */

#include "mex.h"
#include "math.h"
#include "complex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

int f, k, Nloop, Nfreq, kre, m;
double sigma2, gL, C, Delta, VT, VL, Vth, Vlb, dV, Vr, tref, tau_x, Vx, gx, V, gamx, w, pi, byt0;
double u1, r0, mu_in, x0;
double _Complex xibytu,xibytx,x0bytu,x0bytx;
double *params, *xi, *P0, *p0, *freq, *Ar, *Ai, *x1r, *x1i, *V1r, *V1i;

        
/* get inputs */
params = mxGetPr(prhs[0]);
x0 = mxGetScalar(prhs[1]);
mu_in = mxGetScalar(prhs[2]);
sigma2 = mxGetScalar(prhs[3]);
xi = mxGetPr(prhs[4]);
u1 = mxGetScalar(prhs[5]);
r0 = mxGetScalar(prhs[6]);
P0 = mxGetPr(prhs[7]);
p0 = mxGetPr(prhs[8]);
freq = mxGetPr(prhs[9]);
Nfreq = mxGetN(prhs[9]);
m = mxGetM(prhs[9]);

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
tau_x = params[10];
Vx = params[11];
gx = params[12];

Nloop = (int)floor((Vth-Vlb)/dV); /* number of bins */
kre = (int)round((Vr-Vlb)/dV); /* index of reset potential */
pi=3.14159265;
gamx = gx/gL;

/* allocate outputs */
plhs[0] = mxCreateDoubleMatrix(Nfreq,1,mxCOMPLEX);
plhs[1] = mxCreateDoubleMatrix(Nfreq,1,mxCOMPLEX);
plhs[2] = mxCreateDoubleMatrix(Nfreq,1,mxCOMPLEX);

Ar = mxGetPr(plhs[0]);
Ai = mxGetPi(plhs[0]);
x1r = mxGetPr(plhs[1]);
x1i = mxGetPi(plhs[1]);
V1r = mxGetPr(plhs[2]);
V1i = mxGetPi(plhs[2]);

/* initialize arrays for integration */
double *I0,*G,*alpha,*beta;
double _Complex *ru,*rx,*jr,*ju,*jx,*pr,*pu,*px,*psp,*Pu,*Px,*P1;

ru=mxMalloc(Nfreq*sizeof(double _Complex));
rx=mxMalloc(Nfreq*sizeof(double _Complex));
jr=mxMalloc(Nloop*sizeof(double _Complex));
ju=mxMalloc(Nloop*sizeof(double _Complex));
jx=mxMalloc(Nloop*sizeof(double _Complex));
pr=mxMalloc(Nloop*sizeof(double _Complex));
pu=mxMalloc(Nloop*sizeof(double _Complex));
px=mxMalloc(Nloop*sizeof(double _Complex));
psp=mxMalloc(Nloop*sizeof(double _Complex));
Pu=mxMalloc(Nloop*sizeof(double _Complex));
Px=mxMalloc(Nloop*sizeof(double _Complex));
P1=mxMalloc(Nloop*sizeof(double _Complex));
I0=mxMalloc(Nloop*sizeof(double));
G=mxMalloc(Nloop*sizeof(double));
alpha=mxMalloc(Nloop*sizeof(double));
beta=mxMalloc(Nloop*sizeof(double));

V=Vlb;

for (k=0; k<Nloop; k++){
     
    I0[k] = gL*(VL-V)+gx*x0*(Vx-V)+mu_in+gL*Delta*exp((V-VT)/Delta);
    G[k] = -I0[k]/(gL*sigma2);
    alpha[k] = exp(dV*G[k]);
    if (G[k]<.000001) beta[k] = dV/sigma2;
    else beta[k] = (alpha[k]-1)/(G[k]*sigma2);
    
    jr[k] = 0;
    ju[k] = 0;
    jx[k] = 0;
    pr[k] = 0;
    pu[k] = 0;
    px[k] = 0;
    psp[k] = 0;
    Pu[k] = 0;
    Px[k] = 0;
    P1[k] = 0;
    
    V += dV;
    
}


for (f=0; f<Nfreq; f++) {
    /* initialize flux and probability to 0 */
    jr[Nloop-1] = 1+0*I; // boundary condition
    ju[Nloop-1] = 0+0*I;
    jx[Nloop-1] = 0+0*I;
    pr[Nloop-1] = 0+0*I;
    pu[Nloop-1] = 0+0*I;
    px[Nloop-1] = 0+0*I;
    
    w = freq[f]*2*pi;
    
    V = Vth;
    
    for (k = Nloop-2; k>=0; k--){
        
        /* boundary conditions */
        jr[k] = jr[k+1]+dV*I*w*pr[k+1]; if (k==kre) jr[k] -= cexp(-I*w*tref);
        pr[k] = pr[k+1]*alpha[k+1]+C/gL*jr[k+1]*beta[k+1];
        
        /* inhomogenous component - current modulations */
        ju[k] = ju[k+1]+dV*I*w*pu[k+1];
        pu[k] = pu[k+1]*alpha[k+1]+(C/gL*ju[k+1]-u1/gL*(r0*p0[k+1]))*beta[k+1];
        
        /* inhomogenous component - voltage-activated conductance modulations */
        jx[k] = jx[k+1]+dV*I*w*px[k+1];
        px[k] = px[k+1]*alpha[k+1]+(C/gL*jx[k+1]+(V-Vx)*(r0*p0[k+1]))*beta[k+1];
        
        V -= dV;
    }
    
    /* normalize and compute first-order activity */
    
    ru[f] = -ju[0]/jr[0];
    rx[f] = -jx[0]/jr[0];
    
    V = Vlb;
    xibytu = 0; xibytx = 0; x0bytu = 0; x0bytx = 0;
    for (k=0; k<Nloop; k++) {
        
        if (k<kre) psp[k] = 0;
        else psp[k] = cexp(-I*w*tref*((Vth-V)/(Vth-Vr)))/(Vth-Vr);
        
        Pu[k] = pu[k]+ru[f]*(pr[k]+tref*psp[k]);
        Px[k] = px[k]+rx[f]*(pr[k]+tref*psp[k]);
        
        xibytu += (Pu[k]*xi[k]);
        xibytx += (Px[k]*xi[k]);
        x0bytu += Pu[k];
        x0bytx += Px[k];
        byt0 += P0[k];
        
        V += dV;
    }
    
    xibytu = xibytu*dV/tau_x;
    xibytx = xibytx*dV/tau_x;
    x0bytu = x0bytu*x0*dV/tau_x;
    x0bytx = x0bytx*x0*dV/tau_x;
    byt0 = byt0*dV/tau_x;
    
    x1r[f] = creal((xibytu-x0bytu)/(byt0+I*w-gamx*(xibytx-x0bytx)));
    x1i[f] = cimag((xibytu-x0bytu)/(byt0+I*w-gamx*(xibytx-x0bytx)));
    Ar[f] = creal(ru[f]+gamx*(x1r[f]+x1i[f])*rx[f]);
    Ai[f] = cimag(ru[f]+gamx*(x1r[f]+x1i[f])*rx[f]);
    
    
    V = Vlb;
    for (k=0; k<Nloop; k++) {
        P1[k] = Pu[k]+(x1r[f]+x1i[f])*Px[k];
        V1r[f] += creal(P1[k]*V);
        V1i[f] += cimag(P1[k]*V);
        V += dV;
    }
    V1r[f] = V1r[f]*dV;
    V1i[f] = V1i[f]*dV;
    
    
}

mxFree(ru);
mxFree(rx);
mxFree(jr);
mxFree(ju);
mxFree(jx);
mxFree(pr);
mxFree(pu);
mxFree(px);
mxFree(psp);
mxFree(Pu);
mxFree(Px);
mxFree(P1);
mxFree(I0);
mxFree(G);
mxFree(alpha);
mxFree(beta);

}

