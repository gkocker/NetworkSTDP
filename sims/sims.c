#include "mex.h"
#include "fastexponent.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

/* Constants for ran2 random number generator. From Numerical Recipes. */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* declare functions used other than mexFunction*/
double ran2(long *idum);    /* Ran2 uniform random deviate generator from Numerical Recipes */
double gasdev(long *idum);
int myComparamsisonFunction(const void *x, const void *y);

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{

// declare variables
    double *params,*seed,*spk_pad,*xi_c_tmp,*ccg_out,*rates,*V_init,*x_init,*Isyn_init,*Conn_in,*mu_vec,*sig2_vec,*W,*indWpre,*indWpost,*start_pre,*start_post,*Npre,*Npost,*Wout;
    double gL,Delta,VT,VL,C,Vth,Vr,tref,gx,tau_x,Vx,Vxh,Dx;
    double tauE,tauI,p,wmaxEE,wmaxEI,wmaxIE,wmaxII,Fpot,Fdep,tau_pot,tau_dep,Isyn_in,dt_w,nind2;
    double c,xi,xi_c,cc,ci,tstop,trans,dt,t,w,ton,Iapp,pi;
    double ccg_win_rad;  // Window radius for which to calculate cross-correlation functions (ms)
    double ccg_dt;		// Bin size for cross-correlation functions (ms)
    
    long idum;
    int countM,R,r,Ni,Ne,N,n,ne,ni,i,j,k,Tloop,Transloop,count,t_w,Tloopw,indw,nind;
    
    // Get pointers to input arrays
    params = mxGetPr(prhs[0]);
    seed = mxGetPr(prhs[1]);  
    mu_vec = mxGetPr(prhs[2]);
    sig2_vec = mxGetPr(prhs[3]);
    Conn_in = mxGetPr(prhs[4]);
    
    idum = seed[0]<0?(long)seed[0]:-(long)seed[0]; /* Ran2 needs a negative value for a seed. If the seed is not negative, make it negative. */
    
    /* Get paramsameters */
    gL = params[0];
    Delta = params[1];
    VT = params[2];
    VL = params[3];
    C = params[4];
    Vth = params[5];
    Vr = params[6];
    tref = params[7];
    gx = params[8];
    tau_x = params[9];
    Vx = params[10];
    Vxh = params[11];
    Dx = params[12];
    c = params[13];
    tstop = params[14];
    trans = params[15];
    dt = params[16];
    Ni = params[17];
    Ne = params[18];
    tauE = params[19];
    tauI = params[20];
//     p = params[21];
    wmaxEE = params[22];
    wmaxEI = params[23];
    wmaxIE = params[24];
    wmaxII = params[25];
    Fpot = params[26];
    Fdep = params[27];
    tau_dep = params[28];
    tau_pot = params[29];       
    w = params[30];
    ton = params[31];
    ccg_win_rad = params[32];
    ccg_dt = params[33];
    dt_w = params[34];
        
    Tloop = (int)(tstop/dt);
    Tloopw = (int)((tstop-trans)/dt_w);
    Transloop = (int)(trans/dt);
    N = Ne+Ni;
    countM = 50*N*tstop/1000; // max number of spikes to output for spectral calculations (prefactor should be greater than firing rate)
  
//     mexPrintf("%d\n",nind);
    
    /* Set up[ cross-covariance calculations.  Adapted from James Trousdale, Univ. Houston Dept. Mathematics */
    
    /* The following sets of variables are used in the processing of the cross-correlation functions.
	   The gist of the method is this: the spike histories of all cells are retained for the duration
	   of the ccg window radius. After every simulation time step, each cell of each pair in the list
	   of pairs for which we are calculating CCGs is checked for a spike during the preceeding step.
	   If a cell spiked, its partners spike history is checked, and we increment the one-sided
	   un-normalized correlogram corresponding to the past for the pair. By doing this for both cells
	   in each pair, we construct two one-sided correlograms which may be attached to form the full
	   correlogram, which can then be normalized to acquire the cross-correlation function. Additionally,
	   we rotate through the spike history vectors in place to eliminate the need for costly read/write
	   operations. - James*/
    

    
    int num_ccg_pairs;
    num_ccg_pairs = 3;

    int ccg_inds[num_ccg_pairs][2];
    
    /* input neuron indices for pairs to calculate ccgs of */
   ccg_inds[0][0] = 0;
   ccg_inds[0][1] = 1;     // for two-cell network, indices are cells 0 and 1
    
   ccg_inds[1][0] = 0;
   ccg_inds[1][1] = 2;
    
   ccg_inds[2][0] = 1;
   ccg_inds[2][1] = 2;
    
    
	int bins = (int)ceil(ccg_win_rad/ccg_dt); // get ccg_win_rad and ccg_dt from parameters input vector
	int dt_ratio = (int)ceil(ccg_dt/dt);
	int spike_hist_steps = (int)ceil((ccg_win_rad+ccg_dt/2)/dt)+1; // Number of time bins for which it is necessary														   // to retain spike history for the calculation of
																	   // the one-sided correlograms.
	int ind1, ind2;     // neuron indices for calculating cross-correlograms
    
	// spike_hist retains the spike histories of every cell for use in calculation of one-sided correlograms.
    int spike_hist[N][spike_hist_steps];
	for(i = 0; i < N; ++i) {
        for(j = 0; j < spike_hist_steps; ++j) spike_hist[i][j] = 0;
    }
    
    // ccg_o is of size (num_ccg_pairs)x2x(bins+1), and holds the one-sided un-normalized CCGs
    double ccg_o[num_ccg_pairs][2][bins+1];
    for(i = 0; i < num_ccg_pairs; ++i) {
        for(j = 0; j < 2; ++j) {
            for(k = 0; k < bins+1; ++k) ccg_o[i][j][k] = 0;
        }
    }
    
	
    int sh_ind = 0; // sh_ind will track which time bin of the correlogram spike-histories we are currently centered on.
    
    // The ccg_o_inds vectors are essentially "convenience" variables which give the CORRELOGRAM bin corresponding to
	// a time lag given in SIMULATION time bins. In other words, maps a time difference in sim time bins to a time
	// difference in correlogram bins.
	
//     int* ccg_o_inds_1 = new int[spike_hist_steps];
    int ccg_o_inds_1[spike_hist_steps];
    for(i = 0; i < spike_hist_steps; ++i) {
        ccg_o_inds_1[i] = (int)ceil(((double)(i-(double)dt_ratio/2))/dt_ratio);
    }
        
//     int* ccg_o_inds_2 = new int[spike_hist_steps-2];
    int ccg_o_inds_2[spike_hist_steps-2];
    for(i = 0; i < spike_hist_steps-2; ++i) {
        ccg_o_inds_2[i] = (int)ceil(((double)(i-(double)dt_ratio/2+2))/dt_ratio);
    }
    
    /* End code adapted from James */

    double *V, *V0, *tlast, *D_vec, *syn_pre, *syn_pre0, *syn_post, *syn_post0, *Isyn, *Isyn0, *taus, *num_spikes;
    
    
    V = mxMalloc(N*sizeof(double));
    V0 = mxMalloc(N*sizeof(double));
    tlast = mxMalloc(N*sizeof(double));
    D_vec = mxMalloc(N*sizeof(double));
    syn_pre = mxMalloc(N*sizeof(double));
    syn_pre0 = mxMalloc(N*sizeof(double));
    syn_post = mxMalloc(N*sizeof(double));
    syn_post0 = mxMalloc(N*sizeof(double));
    Isyn = mxMalloc(N*sizeof(double));
    Isyn0 = mxMalloc(N*sizeof(double));
    taus = mxMalloc(N*sizeof(double));
    num_spikes = mxMalloc(N*sizeof(double));
    
    double *Conn = mxMalloc(N*N*sizeof(double));
    double *Conn0 = mxMalloc(N*N*sizeof(double));
//     double Conn[N][N];
//     double Conn0[N][N];
    
//     double V[N], V0[N], tlast[N], D_vec[N];
//     double syn_pre[Ne], syn_pre0[Ne], syn_post[Ne], syn_post0[Ne], Isyn[N], Isyn0[N];
//     double Conn[N][N], Conn0[N][N], taus[N];
//     double num_spikes[N];
    
    int indw2[3] = {0,0,0};
    int dims[3] = {N,N,Tloopw};
    
//     mexPrintf("%d\n",Tloopw);
    
    /* Create an mxArray for the output data */
    plhs[0] = mxCreateDoubleMatrix(2,countM,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,Tloop,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(num_ccg_pairs,2*bins+1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(N,1,mxREAL);
    plhs[4] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    
    /* Create a pointer to the output data */
    spk_pad = mxGetPr(plhs[0]);
    xi_c_tmp = mxGetPr(plhs[1]);
    ccg_out = mxGetPr(plhs[2]);
    rates = mxGetPr(plhs[3]);
    Wout = mxGetPr(plhs[4]);
        
    // Get initial conditions
    for(n = 0; n < N; n++){ 
        
        D_vec[n] = gL*sqrt(2*sig2_vec[n]*C/gL);
        V0[n] = gasdev(&idum)*(Vth-Vr);
        tlast[n] = -tstop;
        Isyn[n] = 0;
        Isyn0[n] = 0;
        rates[n] = 0;
        num_spikes[n] = 0;
        syn_pre[n] = 0;
        syn_pre0[n] = 0;
        syn_post[n] = 0;
        syn_post0[n] = 0;
//         
        if(n<Ne) {
            taus[n] = tauE;
        }
        else {
            taus[n] = tauI;
        } 
        
        for (j=0; j<N; j++){
//             Conn0[n][j] = Conn_in[N*j+n];
//             Conn[n][j] = Conn0[n][j];
            Conn0[N*j+n] = Conn_in[N*j+n];
            Conn[N*j+n] = Conn_in[N*j+n];
        }
        
    }
//    
// 
    
    /* Initialize simulation variables */
    t = 0; 
    t_w = 0;
    count = 0;
    cc = sqrt(c);
    ci = sqrt(1-c);
    pi = 3.14159265;
    double sqrtdt = sqrt(dt);
    
//     mexPrintf("%f\n",Fdep);
    
    for (i = 0; i < Tloop; i++) {
        
        t = t+dt;
        xi_c_tmp[i] = gasdev(&idum);       
  
        for (ne = 0; ne < Ne; ne++) { //E cells

            Isyn[ne] = Isyn0[ne]+dt*(-Isyn0[ne]/taus[ne]);
            syn_pre[ne] = syn_pre0[ne]+dt*(-syn_pre0[ne]/tau_pot);
            syn_post[ne] = syn_post0[ne]+dt*(-syn_post0[ne]/tau_dep);

            
            xi_c = cc*D_vec[ne]*xi_c_tmp[i];
            xi = ci*D_vec[ne]*gasdev(&idum);
            V[ne]  =  V0[ne]+(dt/C)*(gL*(VL-V0[ne])+mu_vec[ne]+Isyn[ne]+gL*Delta*EXP((V0[ne]-VT)/Delta))+(sqrtdt/C)*(xi_c+xi);

            if (t-tlast[ne] < tref) {
                V[ne]  =  Vth-(t-tlast[ne])/tref*(Vth-Vr); // linear spike shape from threshold to reset
            }
            
            if(V[ne] >= Vth) {
                
                V[ne] = Vth;                    // ensure V doesn't cross recording threshold (uniform spike shape)
                tlast[ne] = t;                  // update last spike time for refractory period
                syn_pre[ne] = syn_pre[ne]+1;    // update this neuron's presynaptic trace
                syn_post[ne] = syn_post[ne]+1;  // update this neuron's postsynaptic trace
                
                if (t < trans){
                    for (k = 0; k < N; k++){
                        if  (Conn_in[N*ne+k] != 0){
                        Isyn[k] += Conn[N*ne+k];
                        Isyn0[k] = Isyn[k];
//                             Isyn[k] += Conn[k][n];
//                             Isyn0[k] = Isyn[k];
                        }
                    }
                }
                else{
                    
                    num_spikes[ne] += 1;			// Increment the spike count of cell j
                    spike_hist[ne][sh_ind] = 1;		// Set the spike to be counted in to the one-sided correlograms involving this cell.
                    
                    for (k = 0; k < N; k++){
      
                        if (Conn_in[N*ne+k] != 0){
                            if (Conn0[N*ne+k] > 0){
                                Isyn[k] += Conn[N*ne+k];
                                Isyn0[k] = Isyn[k];
//                                 Conn[k][ne] = fmax(Conn0[k][ne]+Fdep*syn_post[k],0);
//                                 Conn0[k][ne] = Conn[k][ne];
                                Conn[N*ne+k] = fmax(Conn0[N*ne+k]+Fdep*syn_post[k],0);
                                Conn0[N*ne+k] = Conn[N*ne+k];
                            
                            }
                        }
                        if (Conn_in[N*k+ne] != 0 && Conn0[N*k+ne] < wmaxEE){
//                             mexPrintf("%d\n",k);
//                             Conn[ne][k] = fmin(Conn0[ne][k]+Fpot*syn_pre[k],wmaxEE);
//                             Conn0[ne][k] = Conn[ne][k];
                                Conn[N*k+ne] = fmin(Conn0[N*k+ne]+Fpot*syn_pre[k],wmaxEE);
                                Conn0[N*k+ne] = Conn[N*k+ne];
                        }
                        
                    }
//                     for (k=(int)start_post[ne]-1; k<(int)start_post[ne+1]-1; k++) {
//                         Isyn[(int)Npost[k]-1] += W[(int)indWpost[k]-1]; // neuron ne's postsynaptic targets get their input currents kicked
//                         W[(int)indWpost[k]-1] += Fdep*syn_post[(int)Npost[k]-1]; // synapses of neuron ne's postsynaptic targets depress
//                         if (W[(int)indWpost[k]-1] < 0)
//                             W[(int)indWpost[k]-1] = 0;
//                     }
//                     for (k=(int)start_pre[ne]-1; k<(int)start_pre[ne+1]-1; k++) {
//                         W[(int)indWpre[k]-1] += Fpot*syn_pre[(int)Npre[k]-1]; // synpases of neuron ne's presynaptic inputs potentiate
//                         if (W[(int)indWpre[k]-1] > wmaxEE)
//                             W[(int)indWpre[k]-1] = wmaxEE;
//                     }
                }
            }
//                /* OLD STDP */
// // // //                 for(j = 0;j<N;j++){ // additive, anti-Hebbian stdp
// // // // //                                         
// // // //                     if (Conn_init[N*ne+j] !=0 && Conn0[j][ne] < wmax[j][ne]){
// // // // // //                     if (Conn_init[j*N+ne] !=0 && Conn0[j][ne] < wmax[j][ne]){
// // // //                         Conn[j][ne] = Conn0[j][ne]+Fdep*syn_post[j]; // synapses neuron ne is presynaptic in potentiate prop. to Fdep(column ne)                    
// // // //                     }
// // // //                     if (Conn_init[N*j+ne] !=0 && Conn0[ne][j] > 0){
// // // // // //                     if (Conn_init[ne*N+j] !=0 && Conn0[ne][j] > 0){
// // // //                         Conn[ne][j] = Conn0[ne][j]+Fpot*syn_pre[j];// synapses neuron ne is postsynaptic in depress (row ne) proportionally Fpot and to their presynaptic trace
// // // //                     }
// // // // //                     
// // // // //                    
// // // //                     if (Conn[j][ne] > wmax[j][ne]){ // check that didn't pass upper bound
// // // // //                         Conn[j][ne] = wmax[j][ne]; 
// // // //                     }
// // // //                     if (Conn[ne][j] < 0){ // check that didn't pass lower bound
// // // // //                         Conn[ne][j] = 0;
// // // //                     }
// // // // //                     
// // // //                     Conn0[j][ne] = Conn[j][ne];
// // // //                     Conn0[ne][j] = Conn[ne][j];
// // // //                 }
// // // //                 Conn0[ne][ne]  =  0;
// // //                                
// // 
// // // //                 for(j = 0;j<N;j++){ // multiplicative stdp
// // // //                     
// // // //                     if (Conn_init[N*ne+j] !=0){
// // // //                         Conn[j][ne] = Conn0[j][ne]+Fdep*Conn0[j][ne]*syn_post[j]; // synapses neuron ne is presynaptic in depress(column ne) proportionally to their postsynaptic trace
// // // //                         Conn0[j][ne] = Conn[j][ne];
// // // //                     }
// // // //                     if (Conn_init[N*j+ne] !=0){
// // // //                         Conn[ne][j] = Conn0[ne][j]+Fpot*(wmax[ne][j]-Conn0[ne][j])*syn_pre[j];// synapses neuron ne is postsynaptic in potentiate (row ne) proportionally to their presynaptic trace
// // // //                         Conn0[ne][j] = Conn[ne][j];
// // // //                     }
// // // //                 }
// // // //                 Conn0[ne][ne]  =  0;
// // //             }
// // //             }
            // update previous-time-step values
            V0[ne] = V[ne];
            Isyn0[ne] = Isyn[ne];
            syn_pre0[ne] = syn_pre[ne];
            syn_post0[ne] = syn_post[ne];
            
            
//             if ( i%(Tloop/Tloopw) == 0){
//                 for (k=(int)start_post[ne]-1; k<(int)start_post[ne+1]-1; k++)
//                     Wout[t_w*nind + (int)(indWpost[k]-1)] = W[(int)indWpost[k]-1];
//             }
//             
        } //end neuron loop 
        
        /* Store Connectivity */
        j = i-Transloop;
        if (j>=0 && j%((Tloop-Transloop)/Tloopw) == 0){
            indw2[2] = t_w;
            for (n=0; n<N; n++){
                    indw2[0] = n;
                for (k=0; k<N; k++){
                    indw2[1] = k;
                    indw = mxCalcSingleSubscript(plhs[4],3,indw2);
//                     Wout[indw] = Conn[n][k];
                    Wout[indw] = Conn[N*k+n];
                }
            }
            
            t_w += 1;
//             mexPrintf("%d\n",j);
//             mexPrintf("%d\n",t_w);
        }
        
        /* Adapted from James Trousdale, Univ. Houston Dept. Mathematics */
        // If cells of interest spiked, add to ccg
        for (j = 0; j < num_ccg_pairs; ++j) {
            ind1 = ccg_inds[j][0];
            ind2 = ccg_inds[j][1];
            
            // Check the first cell in pair j. If it spiked, then search over the recent spike history of its partner. If its partner spiked in a bin a
            // certain number of simulation bins in the past, then increment the appropriate correlogram bin.
            if(spike_hist[ind1][sh_ind] == 1)
                for( k = 0; k < spike_hist_steps; ++k)
                    if(spike_hist[ind2][(sh_ind+k)%spike_hist_steps] == 1)
                        ccg_o[j][0][ccg_o_inds_1[k]] += 1;
            
            // Check the second cell in pair j. If it spiked, then search over the recent spike history of its partner. If its partner spiked in a bin a
            // certain number of simulation bins in the past, then increment the appropriate correlogram bin. Note that we avoid double counting spikes
            // occuring simultaneously (same simulation time bin).
            if(spike_hist[ind2][sh_ind] == 1)
                for( k = 1; k < spike_hist_steps-1; ++k)
                    if(spike_hist[ind1][(sh_ind+k)%spike_hist_steps] == 1)
                        ccg_o[j][1][ccg_o_inds_2[k-1]] += 1; 
        }
         
        for (j = 0; j < N; ++j) {
            spike_hist[j][(sh_ind-1+spike_hist_steps)%spike_hist_steps] = 0; // Erase spikes from the correlogram history which move outside of the CCG window radius.
        }
        sh_ind = (sh_ind-1+spike_hist_steps)%spike_hist_steps; // Incrememnt the spike history bin we are centered on.
        /* End code adapted from James */
        
    } //end time loop
    
    
//     mexPrintf("%d\n",(int)sizeof(Isyn));
//     for(n=0; n<N; n++)
//     mexPrintf("%f\n",Isyn[n]);
    
    
// calculate rates
// double rates[N];
for (i = 0; i<N; i++)
    rates[i] = (double)num_spikes[i]/tstop; // sp/ms

// calculate two-sided cross-correlograms
double ccg[num_ccg_pairs][2*bins+1];
for (i = 0; i < num_ccg_pairs; ++i) {
    ccg[i][bins] = ccg_o[i][0][0]+ccg_o[i][1][0]; // 0 lag ccg
    for(j = 0; j < bins; ++j) {
        ccg[i][bins-(j+1)] = ccg_o[i][1][j+1];
        ccg[i][bins+(j+1)] = ccg_o[i][0][j+1];
    }
}

// normalize into cross-covariance functions
for (i = 0; i < num_ccg_pairs; ++i){
    for (j = 0; j < 2*bins+1; ++j){
        ccg[i][j]=rates[ccg_inds[i][1]]*(ccg[i][j]/(num_spikes[ccg_inds[i][1]]*ccg_dt)-rates[ccg_inds[i][0]]);
    }
}

// copy into output matrix
for (i = 0; i < num_ccg_pairs; ++i){ 
    for (j = 0; j < 2*bins+1; ++j){
//             ccg_out[(2*bins+1)*i+j] = ccg[i][j];        
        ccg_out[num_ccg_pairs*j+i] = ccg[i][j];
    }
}
        

mxFree(V);
mxFree(V0);
mxFree(tlast);
mxFree(D_vec);
mxFree(syn_pre);
mxFree(syn_pre0);
mxFree(syn_post);
mxFree(syn_post0);
mxFree(Isyn);
mxFree(Isyn0);
mxFree(taus);
mxFree(num_spikes);
mxFree(Conn);
mxFree(Conn0);


} //end mex function
      
int myComparamsisonFunction(const void *x, const void *y) {

    // x and y are pointers to doubles.

    // Returns -1 if x < y
    //          0 if x  ==  y
    //         +1 if x > y

    double dx, dy;

    dx  =  *(double *)x;
    dy  =  *(double *)y;

    if (dx < dy) {
        return -1;
    } 
    else if (dx > dy) {
        return +1;
    }
    return 0;
}


double ran2(long *idum)
{
 int j;
 long k;
 static long idum2=123456789;
 static long iy=0;
 static long iv[NTAB];
 float temp;

 if (*idum <= 0) {
  if (-(*idum) < 1) *idum=1;
  else *idum = -(*idum);
  idum2=(*idum);
  for (j=NTAB+7;j>=0;j--) {
   k=(*idum)/IQ1;
   *idum=IA1*(*idum-k*IQ1)-k*IR1;
   if (*idum < 0) *idum += IM1;
   if (j < NTAB) iv[j] = *idum;
  }
  iy=iv[0];
 }
 k=(*idum)/IQ1;
 *idum=IA1*(*idum-k*IQ1)-k*IR1;
 if (*idum < 0) *idum += IM1;
 k=idum2/IQ2;
 idum2=IA2*(idum2-k*IQ2)-k*IR2;
 if (idum2 < 0) idum2 += IM2;
 j=iy/NDIV;
 iy=iv[j]-idum2;
 iv[j] = *idum;
 if (iy < 1) iy += IMM1;
 if ((temp=AM*iy) > RNMX) return RNMX;
 else return temp;
}

double gasdev(long *idum)
{
    double ran2(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;

    if (*idum < 0) iset=0;
    if  (iset == 0) {
        do {
            v1=2.0*ran2(idum)-1.0;
            v2=2.0*ran2(idum)-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    } else {
        iset=0;
        return gset;
    }
}