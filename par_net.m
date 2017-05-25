%%%Network parameters
NE=1000; 
NI=0;
N=NE+NI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate adjacency matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pEE=.15;
pref = [pEE*ones(N,N)];
adj = rand(N); adj(adj>pref) = 0; adj(adj~=0) = 1;

%%% embed three neuron microcircuit with direct, reciprocal and common
%%% input motifs as in Fig. 2
adj(1:3,1:3) = 0;
adj(1,2) = 1; adj(2,1) = 1; adj(3,1) = 1;

adj(1:N+1:end) = 0; % remove autapses

%%% calculate motif frequencies
p0 = sum(sum(adj))/(N^2);
q_rec0 = sum(sum(adj.*adj.'))/N^2 - p0^2;
q_div0 = sum(sum(adj*adj.'))/(NE^3) - p0^2;
q_con0 = sum(sum(adj.'*adj))/(NE^3) - p0^2;
q_ch0 = sum(sum(adj*adj))/(NE^3) - p0^2;

wmaxE = 5/(N*p0);
wmaxI = 0;
wmax = zeros(N); wmax(:,1:NE) = wmaxE; wmax(:,NE+1:N) = wmaxI;

p_cond = wmaxE*.25;
Conn(1:NE,1:NE)=p_cond*adj(1:NE,1:NE); % to E from E

%%% asymmetric weights for the three-neuron microcircuit (as in Fig. 3)
Conn(1,2) = wmaxE*.3;
Conn(2,1) = wmaxE*.7;
Conn(3,1) = wmaxE*.3;

%%% calculate initial motif strengths
p_init = sum(sum(Conn))/(NE^2);
q_div_init = sum(sum(Conn*Conn.'))/(NE^3)-p_init^2;
q_con_init = sum(sum(Conn.'*Conn))/(NE^3)-p_init^2;
q_ch_init = sum(sum(Conn*Conn))/(NE^3)-p_init^2;
q_rec_init = sum(sum(Conn.*Conn.'))/(NE^2)-p_init^2;
q_Xrec_init = sum(sum(adj.*Conn.'))/(NE^2)-p_init*p0;
q_Xdiv_init = sum(sum(Conn*adj.'))/(NE^3)-p_init*p0;
q_Xcon_init = sum(sum(Conn.'*adj))/(NE^3)-p_init*p0;
q_XchA_init = sum(sum(Conn*adj))/(NE^3)-p_init*p0;
q_XchB_init = sum(sum(adj*Conn))/(NE^3)-p_init*p0;

%%% synaptic time constants
%%% note: sims are not set up for synaptic delays yet
taud=0*ones(N,1); %synaptic delays for output of each cell, ms
tauE=5; tauI=2;
taus=zeros(N,1);
taus(1:NE)=tauE; taus(NE+1:N)=tauI; % time constants for synaptic output of each cell, ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% specify the stdp rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shift_stdp = 0; %rightward shift is positive, leftward shift is negative

%%% Hebbian and slow pot/fast dep anti-Hebbian
tau_pos = 15; tau_neg = tau_pos*2; %roughly consisent with Froemke & Dan 2002 (Hebbian)

deleps = -wmaxE/50000; %for temporally asymmetric
Fneg = -wmaxE/5000;
Fpos = (-Fneg*tau_neg+deleps)/tau_pos;

tau_stdp  =  zeros(N); 
tau_stdp(~logical(tril(tau_stdp)))  =  tau_neg; %post
tau_stdp(~logical(triu(tau_stdp)))  =  tau_pos; %pre