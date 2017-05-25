%%% Euler method for self-consistent changes in synaptic
%%% weights and spiking cross-covariance functions
%%% first call par_net.m and par_cell.m to define parameters

tic;

%%% compile C code for calculating marginal statistics
mex theory_EIFresp.c
mex theory_EIFfpt.c
mex thin_x.c

% dt_w = 100000; % long timescale time step
% Tmax_w = 10000000; % total time for plasticity
dt_w = 1;
Tmax_w = 1;

c_in = 0; % correlation of external input

num_loop_stdp = ceil(Tmax_w/dt_w)+1; 

Conn_t = cell(num_loop_stdp,1); % to store connectivity
Conn_t{1} = sparse(Conn); % initial condition for weights

p = zeros(num_loop_stdp,1); % to store motif statistics
q_div = zeros(num_loop_stdp,1);
q_con = zeros(num_loop_stdp,1);
q_ch = zeros(num_loop_stdp,1);

p(1) = sum(sum(Conn))/(N^2);
q_div(1) = sum(sum(Conn*Conn.'))/(N^3)-p(1)^2;
q_con(1) = sum(sum(Conn.'*Conn))/(N^3)-p(1)^2;
q_ch(1) = sum(sum(Conn*Conn))/(N^3)-p(1)^2;

params = zeros(26);
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
params(14) = Fpos;
params(15) = Fneg;
params(16) = wmaxE;
params(17) = wmaxI;
params(18) = tau_pos;
params(19) = tau_neg;
params(20) = N;
params(21) = NE;
params(22) = NI;
params(23) = tauE;
params(24) = tauI;
params(25) = shift_stdp;
params(26) = dt_w;


%%% run Euler method
for n = 2:num_loop_stdp
   
    n
    [rates,t_th,yy,dConndt,~,~] = LR_caller(Conn,mu_vec,sig2_vec,params,xi,taus,taud,adj,c_in); % compute cross-cov fxn given weight matrix and change in weight matrix given STDP rule and cross-cov fxn
        
    Conn = Conn+dConndt.*dt_w; % update connectivity
    Conn_t{n} = sparse(Conn); % store connectivity
    
    % calculate motif strengths
    p(n) = sum(sum(Conn))/(N^2); %connectivity / mean synaptic weight of random graph with same number of connections as the network
    q_div(n) = sum(sum(Conn*Conn.'))/(N^3)-p(n)^2;
    q_con(n) = sum(sum(Conn.'*Conn))/(N^3)-p(n)^2;
    q_ch(n) = sum(sum(Conn*Conn))/(N^3)-p(n)^2;   
    
    mean(rates)*1000
    
end

toc;