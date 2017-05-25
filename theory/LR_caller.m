function [rates,t_th,yy,dConndt,freq,yyt] = LR_caller(Conn,mu_vec,sig2_vec,params_in,xi,taus,taud,adj,dE)

%%% Compute cross-spectra and cross-covariance functions for pop of EIF
%%% neurons with given weight matrix
%%% All for renewal EIF models; make sure gx=0 in par_cell (auto-cov. function is only obtainable for renewal spiking)
%%% For details of this method, see Trousdale et al., PLoS Comp. Biol. 2012 

gL = params_in(1);
C = params_in(2);
Delta = params_in(3);
VT = params_in(4);
VL = params_in(5);
Vth = params_in(6);
Vlb = params_in(7);
dV = params_in(8);
Vr = params_in(9);
tref = params_in(10);
tau_x = params_in(11);
Vx = params_in(12);
gx = params_in(13);
Fpos = params_in(14);
Fneg = params_in(15);
wmaxE = params_in(16);
wmaxI = params_in(17);
tau_pos = params_in(18);
tau_neg = params_in(19);
N = params_in(20);
NE = params_in(21);
NI = params_in(22);
tauE = params_in(23);
tauI = params_in(24);
shift_stdp = params_in(25);
dt_w = params_in(26);

D_vec = gL*sqrt(2*sig2_vec*C/gL); % actual noise intensity
wmax = zeros(N); wmax(:,1:NE) = wmaxE; wmax(:,NE+1:N) = wmaxI;

Tmax = 300; % Maximum time lag over which to calculate cross-correlations (ms)
dt = 1; % Bin size for which to calculate cross-correlations (ms)

% Generate a vector of frequencies at which to solve for the spectral
% statistics in order to generate cross-correlations with maximum lag Tmax
% and bin size dt.
df = 1/2/Tmax;
fmax = 1/2/dt;
freq = -fmax:df:(fmax-df);
ind0 = find(abs(freq) < 1e-6);
freq(ind0) = 1e-6;
nfreq = length(freq);

% solve for firing rates via fixed point iteration
params = zeros(13,1);
params(1) = gL;
params(2) = C;
params(3) = Delta;
params(4) = VT;
params(5) = VL;
params(6) = Vth;
params(7) = Vlb;
params(8) = dV;
params(9) = Vr;
params(10) = tref;
params(11) = tau_x;
params(12) = Vx;
params(13) = gx;


mu = mu_vec(1); sigma2 = sig2_vec(1); D = D_vec(1);
[~,~,~,r0,~] = theory0(mu,sigma2,params,xi); %initial estimate
rates = zeros(N,1)+r0; %Hz/1000
rates_temp = zeros(N,1);

%%% multiply weight matrix by net synaptic current for exponential synaptic kernel
syn_int_E = (taus(1:NE)*ones(1,N))'; % ms
Conn = syn_int_E.*Conn; % uA/cm^2 * ms

%%% compute firing rates
num_rate_fp_its = 20;
for i = 1:num_rate_fp_its
        
    for j = 1:N
        mu = mu_vec(j); sigma2 = sig2_vec(j); D = D_vec(j);
        
        mu = mu_vec(j); 
        [~,~,~,rates_temp(j),~] = theory0(mu+Conn(j,:)*rates,sigma2,params,xi);
        
    end
    rates = rates_temp;
end

%%% calculate uncoupled transfer function and Fourier transform of synaptic
%%% kernel

u1 = 1; 
At = zeros(N,nfreq);	 % susceptibility / linear response, sp/s / uA/cm^2
Ft = zeros(N,nfreq);	 % synaptic kernel (Fourier domain)
Ct0 = zeros(N,nfreq);

for i=1:N
    
    Ft(i,:) = exp(1i*-2*pi*freq*taud(i))./(1+1i*2*pi*freq*taus(i));
   
    mu = mu_vec(i); sigma2 = sig2_vec(i); D = D_vec(i);
    mu_eff = mu+Conn(i,:)*rates; % calculate mean of input   
    [At(i,:),~,~,~,~,~,Ct0(i,:)] = theory1(u1,mu_eff,sigma2,params,freq,xi);
   
end

%%% solve matrix equations for cross-spectra of each pair, store in yy
d=0; %input correlation c=d*(u1^2)
yyt = zeros(N,N,nfreq);
I = eye(N);

for j = 1:nfreq
    K = zeros(N,N);
    yy0 = zeros(N,N);
    for k = 1:N
        for l = 1:N
            
            if k <= NE && l <= NE % input correlation only to E cells
                d = dE;
            end

            Pss = d*2*(C/gL/1000)*((u1^2)*D_vec(k)*D_vec(l)); % white common noise

            K(k,l) = Conn(k,l)*At(k,j)*Ft(l,j);
            yy0(k,l) = Pss*At(k,j)*conj(At(l,j));
            
        end

        if k <= NE
            d=dE;
        else d=0;
        end

        Pss = d*2*(C/gL/1000)*((u1^2)*D_vec(k)^2); %white common noise

        yy0(k,k) = Ct0(k,j)+Pss*abs(At(k,j))^2;
    end 
    
    % compute cross-spectra
    yyt(:,:,j) = (I-K)\(yy0)/(I-K');
end


%%% Calculate the inverse Fourier transform of the auto-/cross-spectra for
%%% every E-E pair and store the results backe in yy.

yy = zeros(N,N,nfreq);

for i = 1:N
    for j = 1:N
        [t_th,temp_ccg] = inv_f_trans_on_vector(freq,squeeze(yyt(i,j,:)));
        yy(i,j,:) = real(temp_ccg);
    end
end

%%% compute dConn/dt for E-E connections from convolution of stdp windows and cross-covariances

% transform back to regular weight matrix
Conn = Conn./syn_int_E;

k_0t = find(abs(t_th-shift_stdp)==min(abs(t_th-shift_stdp))); % t_th is vector of lags
dConndt = zeros(N,N); % change in trial-averaged synaptic weight per ms

% STDP RULE: specify for negative time lags
stdp_dep = Fneg.*exp(-abs(t_th(2:k_0t))./tau_neg).'; %additive depression, hard lower bound

% STDP RULE: specify for positive time lags
stdp_pot = Fpos.*exp(-abs(t_th(k_0t:end))./tau_pos).'; %additive potentiation, hard upper bound

for i=1:NE %only plasticity of E-E synapses
    for j=1:NE
        
        if adj(i,j)==1 % only synapses that exist are plastic
            yy_stdp = squeeze(yy(i,j,:)) + rates(i)*rates(j); % spike train cross-correlation
            
            % additive STDP
            dConndt(i,j) = (Conn(i,j)>0)*trapz(stdp_dep.*yy_stdp(2:k_0t))*dt + (Conn(i,j)<wmax(i,j))*trapz(stdp_pot.*yy_stdp(k_0t:end))*dt;

            if abs(Conn(i,j)+dConndt(i,j)*dt_w)>abs(wmax(i,j))
                dConndt(i,j) = (wmax(i,j)-Conn(i,j))/dt_w;
            elseif Conn(i,j)+dConndt(i,j)*dt_w<0 
                dConndt(i,j) = (0-Conn(i,j))/dt_w;
            end

            if Conn(i,j) == wmax(i,j) || Conn(i,j) == 0
                dConndt(i,j) = 0;
            end

        end
    end
end

end %function