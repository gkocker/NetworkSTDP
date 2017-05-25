function [r1 r0 x1 x0 P0 V1 C0] = theory1(u1,mu_in,sigma2,params,freq,xi)

[P0,p0,~,r0,x0] = theory0(mu_in,sigma2,params,xi);
[r1,x1,V1] = theory_EIFresp(params,x0,mu_in,sigma2,xi,u1,r0,P0,p0,freq);
[f0] = theory_EIFfpt(params,mu_in,sigma2,freq);
C0 = r0.*(1+2.*real(f0./(1-f0))); %renewal power spectrum - false

end