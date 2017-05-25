%%% Neuron parameters
Vlb = -100;
Vth = 30; % recording threshold
dV = .1;
V = Vlb+dV:dV:Vth;

SA = 2.5*10^-4;
gL = .1; %mS/sq cm
Delta = 1.4; %originally 1.4
VT = -48; %mV
Vr = -72;
VL = Vr; % mV
C = 1; %uF/sq cm^2
tref = 2;

sig2_vec = zeros(N,1);
sig2_vec(:) = 81;

mu_vec = zeros(N,1);
mu_vec(:) = 1;

%%% parameters for a slow voltage-activated conductance, as in Ocker & Doiron, J. Neurophysiol. 2014
%%% set gx=0 to remove this conductance, which is necessary to calculate the baseline covariance C0 for the linear response theory
%%% the power spectra are computed from renewal theory; adding slow conductances makes the spike train non-renewal

gx = 0; %mS/sq cm, .3 for with
Dx = 8; %slope factor
Vx = -85; %reversal potential of KCNQ
Vxh = -40; %half-activation of KCNQ, -45 / -37
xi = 1./(1+exp(-(V-Vxh)/Dx)); %steady-state activation x_infinity
tau_x = 200; %ms, only measured for control