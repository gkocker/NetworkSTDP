%%% initialization
% clear all
% par

%mex -I/usr/local/include -lgsl sims2.c;

%mex -I/usr/local/include -lgsl -g CFLAGS="\$CFLAGS -std=c99" sims2.c;
mex -g CFLAGS="\$CFLAGS -std=c99" sims.c;

dt=.1; %ms
dt_ds=.5; %ms for spk train statistics
% tstop=5000000;%ms
%trans=100;
%tstop=2e7+trans;

trans = 500;
%tstop = 15000000+trans;
tstop = 5500+trans;
% tstop = 2.2e6;
dt_w=(tstop-trans)/200; %ms for saving Connectivity matrix
if rem((tstop-trans),dt_w) ~= 0
   error('Save points do not divide sim time'); 
end
if tstop/dt_w > tstop/dt
    error('Trying to save too many time points');
end


fs=1000/dt_ds; %Hz, sampling frequency
N_net = 1; % number of adjacency matrices
R=1; % realizations per adjacency matrix

tstop
R

if matlabpool('size')==0 && R>1
	matlabpool open; %initialize cluster
end

c=0
w=0;
ton=tstop+10;

% xcov function
lag = 100; %ms, max lag   for spk train corr fxn, multiple of dt_ds
maxlag = ceil(lag/dt_ds);
xcov_spk = zeros(2*maxlag+1,3); % second dimension: number of pairs
rates = zeros(N,1);

Conn_st = cell(R,N_net);
if N == 2
	Nst = 2;
else
	Nst = 200;
end

p = zeros(N_net,R,ceil((tstop-trans)/dt_w));
p_mon = zeros(N_net,R,ceil((tstop-trans)/dt_w));
p_di = zeros(N_net,R,ceil((tstop-trans)/dt_w));
q_div = zeros(N_net,R,ceil((tstop-trans)/dt_w));
q_con = zeros(N_net,R,ceil((tstop-trans)/dt_w));
q_ch = zeros(N_net,R,ceil((tstop-trans)/dt_w));
q_Xdiv = zeros(N_net,R,ceil((tstop-trans)/dt_w));
q_Xcon = zeros(N_net,R,ceil((tstop-trans)/dt_w));
q_XchA = zeros(N_net,R,ceil((tstop-trans)/dt_w));
q_XchB = zeros(N_net,R,ceil((tstop-trans)/dt_w));

params_st = cell(N_net,1);

tic 

for net = 1:N_net
    
    net
    R_count = R;
    %     par_net;
    % params file for mex
    params = zeros(35,1);
    params(1) = gL;
    params(2) = Delta;
    params(3) = VT;
    params(4) = VL;
    params(5) = C;
    params(6) = Vth;
    params(7) = Vr;
    params(8) = tref;
    params(9) = gx;
    params(10) = tau_x;
    params(11) = Vx;
    params(12) = Vxh;
    params(13) = Dx;
    params(14) = c;
    params(15) = tstop;
    params(16) = trans;
    params(17) = dt;
    params(18) = NI;
    params(19) = NE;
    params(20) = tauE;
    params(21) = tauI;
    params(22) = pEE;
    params(23) = wmaxE; %EE
    params(24) = wmaxI; %EI
    params(25) = wmaxE; %IE
    params(26) = wmaxI; %II
    params(27) = Fpos;
    params(28) = Fneg;
    params(29) = tau_neg;
    params(30) = tau_pos;
    params(31) = w;
    params(32) = ton;
    params(33) = lag;
    params(34) = dt_ds;
    params(35) = dt_w;
    
    params_st{net} = params;
    r_count = 0;

   parfor r=1:R %realizations
        
        % disp(r)
        tic

        %%% run sims
        seed = -round(10000*rand);
        [~,xi_out,xcov_tmp,rates_tmp,Conn_out] = sims(params,seed,mu_vec,sig2_vec,Conn);
        
        xcov_spk = xcov_spk+xcov_tmp';
        rates = rates+rates_tmp;
                
    	Conn_st{r,net} = Conn_out(1:Nst,1:Nst,:); % store subsample of connectivity
        
        %%% calculate motif statistics     
        p_tmp = zeros(1,ceil((tstop-trans)/dt_w));
        p_mon_tmp = zeros(1,ceil((tstop-trans)/dt_w));
        p_di_tmp = zeros(1,ceil((tstop-trans)/dt_w));
        q_div_tmp = zeros(1,ceil((tstop-trans)/dt_w));
        q_con_tmp = zeros(1,ceil((tstop-trans)/dt_w));
        q_ch_tmp = zeros(1,ceil((tstop-trans)/dt_w));
        q_Xrec_tmp = zeros(1,ceil((tstop-trans)/dt_w));
        q_Xdiv_tmp = zeros(1,ceil((tstop-trans)/dt_w));
        q_Xcon_tmp = zeros(1,ceil((tstop-trans)/dt_w));
        q_XchA_tmp = zeros(1,ceil((tstop-trans)/dt_w));
        q_XchB_tmp = zeros(1,ceil((tstop-trans)/dt_w));


        for s=1:floor((tstop-trans)/dt_w)
            
            Conn_tmp = squeeze(Conn_out(:,:,s));
            p_tmp(s) = sum(sum(Conn_tmp))/(N^2);
            q_div_tmp(s) = sum(sum(Conn_tmp*Conn_tmp.'))/(NE^3)-p_tmp(s)^2;
            q_con_tmp(s) = sum(sum(Conn_tmp.'*Conn_tmp))/(NE^3)-p_tmp(s)^2;
            q_ch_tmp(s) = sum(sum(Conn_tmp*Conn_tmp))/(NE^3)-p_tmp(s)^2;
            q_Xrec_tmp(s) = sum(sum(adj.*Conn_tmp.'))/(NE^2)-p_tmp(s)*p0;
            q_Xdiv_tmp(s) = sum(sum(Conn_tmp*adj.'))/(NE^3)-p_tmp(s)*p0;
            q_Xcon_tmp(s) = sum(sum(Conn_tmp.'*adj))/(NE^3)-p_tmp(s)*p0;
            q_XchA_tmp(s) = sum(sum(Conn_tmp*adj))/(NE^3)-p_tmp(s)*p0;
            q_XchB_tmp(s) = sum(sum(adj*Conn_tmp))/(NE^3)-p_tmp(s)*p0;

            Conn_tmp = [];
            
            %             adj = zeros(N);
            %             adj(Conn_out(:,:,s)~=0) = 1;
            %             adj2 = adj - adj.';
            %             p_mon_tmp(s) = sum(sum(Conn_out(:,:,s)(adj2==1)))/(NE^2);
            %
            %             adj2 = adj.*adj.';
            %             p_di_tmp(s) = sum(sum(Conn_out(:,:,s)(adj2==1)))/(NE^2);
            %
        end
        
        p(net,r,:) = p_tmp;
        p_mon(net,r,:) = p_mon_tmp;
        p_di(net,r,:) = p_di_tmp;
        q_div(net,r,:) = q_div_tmp;
        q_con(net,r,:) = q_con_tmp;
        q_ch(net,r,:) = q_ch_tmp;
        q_Xrec(net,r,:) = q_Xrec_tmp;
        q_Xdiv(net,r,:) = q_Xdiv_tmp;
        q_Xcon(net,r,:) = q_Xcon_tmp;
        q_XchA(net,r,:) = q_XchA_tmp;
        q_XchB(net,r,:) = q_XchB_tmp;

        toc

    end % realizations
    
    xcov_spk = xcov_spk./R;
    rates = rates./R;
    
end % networks

toc

if matlabpool('size') ~= 0
    matlabpool close
end
    
% clear spk_pad, clear ind1, clear ind2, clear spk1_tmp, clear spk2_tmp
% if c == 0
% save(strcat(datestr(now,1),'2cell','c0.mat'),'rates','Conn_st','params_st','tstop','trans','dt','dt_w','R','N_net','c');
% else
% save(strcat(datestr(now,1),'2cell','c',num2str(c),'.mat'),'rates','Conn_st','params_st','tstop','trans','dt','dt_w','R','N_net','c');
% end

% if c == 0
% save(strcat('sims_',datestr(now,1),'_2cell','_c0.mat'));
% else
% save(strcat('sims_',datestr(now,1),'_2cell','_c',num2str(c),'.mat'));
% end


%save(strcat('sims_',datestr(now,1),'_xcov1.mat'),'rates','xcov_spk','NE','NI','tstop','trans','dt','R','dt_ds','lag','params','mu_vec','sig2_vec','Conn','params_st');
% save(strcat('sims_NE',num2str(NE),'-p0',num2str(p0),'-p_cond',num2str(p_cond/wmaxE),'.mat'),'p','p_mon','p_di','q_div','q_con','q_ch','tstop','dt_w','rates','xcov_spk',w'p0','p_cond','R','N_net','deleps','params_st');

