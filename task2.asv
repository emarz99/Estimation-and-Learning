%% TASK 2

clear
close all
clc



%% CHOOSE PRBS/SINE-SWEEP


opt_flag = 2;   % SINE-SWEEP

%% PARALLEL POOL
model = 'Simulator_Single_Axis_mod1';
pp = gcp;
if pp.Connected ~= 1
    myPool = parpool;
    myCluster = parcluster('Processes');
    delete(myCluster.Jobs)
end

%% LOAD PARAMETERS (in base workspace)
% Parameters which shall be passed to simulink
parameters_controller
% Initial model (target) parameters
Xu = -0.1067;
Xq =  0.1191;
Mu = -5.9856;
Mq = -2.6463;
Xd = -10.1647;
Md =  450.9085;
theta = [Xu,Xq,Mu,Mq,Xd,Md];
[A,B,C,D] = quad2grey(Xu,Xq,Mu,Mq,Xd,Md); % pass them as Single Input if preparing simulation
% Noise
% noise.Enabler = 0;
noise.Enabler = 1;
noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]
noise.vel_stand_dev = noise.Enabler * 0.01;                                 %[m/s]
noise.attitude_stand_dev = noise.Enabler * deg2rad(0.0076);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(0.01);                   %[rad/s]
% Delays
delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;
%
SetPoint=[0,0];
% sample_time = 0.004;

%% LOAD PARAMETERS

open_system(model);
load_system(model);
MWS = get_param(model,'ModelWorkspace')
 MWS.DataSource = 'MATLAB File';
 MWS.FileName = 'DATA';
 reload(MWS)

%% Setting bounds
lb   = [0.01,11,50];
ub   = [10,20,70];
nvar = 3;


%% GENETIC ALGORITHM
Xu = -0.1067;
Xq =  0.1191;
Mu = -5.9856;
Mq = -2.6463;
Xd = -10.1647;
Md =  450.9085;
g  = -9.81;

UNC_Xu = sqrt(1.4450e-05)*3;
UNC_Xq = sqrt(1.6332e-06)*3;
UNC_Mu = sqrt(0.0037)*3;
UNC_Mq = sqrt(2.5607e-04)*3;
UNC_Xd = sqrt(.5071e-09)*3;
UNC_Md = sqrt(0.0057)*3;

ParaNom = [Xu Xq Mu Mq Xd Md];
ParaStd = [UNC_Xu  UNC_Xq  UNC_Mu UNC_Mq UNC_Xd UNC_Md]/100;

N = 15;
tic;
for i = 1:N
 theta = ParaNom .* (1 + randn.*ParaStd);
 VECT(:,i) = theta;
open_system(model);



load_system(model);
MWS = get_param(model,'ModelWorkspace');
MWS.DataSource = 'MATLAB File';
MWS.FileName = 'DATA';
 reload(MWS)
 
opt1 = optimoptions('ga','PlotFcn','gaplotbestf','UseParallel',false,...
    'UseVectorized',false,'MaxGenerations',8,'PopulationSize',8);
[opt_signal,covScore,exitflag,output,population,scores] = ...
ga(@(ExcitationM) build_input1(ExcitationM,theta),nvar,[],[],[],[],lb,ub,[],opt1);
eta(:,i) = opt_signal;

J(i) = covScore;

grid on
axis padded
fprintf('Genetic algorithm time: %.2fs\n',toc)

end
toc;

%%
N = 15;
MAT_INC = zeros(N,N);
ts = 0.004;
 for i = 1:N
    
     Xu_MC = VECT(1,i);
     Xq_MC = VECT(2,i);
     Mu_MC = VECT(3,i);
     Mq_MC = VECT(4,i);
     Xd_MC = VECT(5,i);
     Md_MC = VECT(6,i);
     
     A     = [Xu_MC, Xq_MC, -9.81;
              Mu_MC, Mq_MC,     0;
              0,     1,         0];

     B      = [Xd_MC; Md_MC;     0];

     C      = [1,     0,     0;
              0,     1,     0;
              0,     0,     1;
              Xu_MC, Xq_MC, 0]; 

     D      = [0;0;0;Xd_MC];

     tic;
     for j = 1:N
         f1               = eta(1,j);
         f2               = eta(2,j);
         T                = eta(3,j);
         simulation_time  = T;
         d_T              = (0:ts:simulation_time);
         NN               = length(d_T);
         ExcitationM      = zeros(NN,2);
         ExcitationM(:,1) = d_T;
         
         f                = f1+(f2-f1)/T.*d_T;  
         u_t              = 0.1 * sin( 2 *pi.*f.*d_T);
         ExcitationM(:,2) = u_t;
         simOut           = sim ('Simulator_Single_Axis1');

         % Input Controllato
         Mt = simOut.Mtot(1:end);

         % Output
         ax   = simOut.ax(1:end);
         q    = simOut.q(1:end);
    
      % Definizione sistema idgrey
parametres = {'Xu',Xu_MC;'Xq',Xq_MC;'Mu',Mu_MC;'Mq',Mq_MC;'Xd',Xd_MC;'Md',Md_MC;};

odefun = @quad2grey;
fcn_type = 'c';
sys0 = idgrey(odefun,parametres,fcn_type);   
         % Dataset Definition
       data = iddata( [q,ax],Mt,ts);
data.InputName = 'Input';
data.OutputName = {'q','ax'};

%trasformazione Dataset da time-domain to frequency domain
DATA_f = fft(data);

% Grey-box estimation
OPT = greyestOptions('SearchMethod','auto');
sys = greyest(DATA_f,sys0,OPT);

UNC_M  =  getcov(sys);
tr     =  trace(UNC_M);

 MAT_INC(i,j) = tr;
         
     end
     toc;
     
 end
 
 %% Calcolo Media e minimo

  MEDIA_MIN = 10;
 MAX_MIN   = 10;
 for k = 1: N
     MEDIA  (k) = mean(MAT_INC(:,k));
     MASSIMO(k) = max (MAT_INC(:,k));
     
     if (MEDIA(k) < MEDIA_MIN)
         MEDIA_MIN = MEDIA(k);
         k_min     = k;
     end
     if(MASSIMO(k)< MAX_MIN)
         MAX_MIN = MASSIMO(k);
         k_max = k;
     end
 end
%%
Xu = -0.1067;
Xq =  0.1191;
Mu = -5.9856;
Mq = -2.6463;
Xd = -10.1647;
Md =  450.9085;
g  = -9.81;

 A     = [Xu, Xq, -9.81;
              Mu, Mq,     0;
              0,     1,         0];

     B      = [Xd; Md;     0];

     C      = [1,     0,     0;
              0,     1,     0;
              0,     0,     1;
              Xu, Xq, 0]; 

     D      = [0;0;0;Xd];

      f1               = eta(1,6);
         f2               = eta(2,6);
         T                = eta(3,6);
         simulation_time  = T;
         d_T              = (0:ts:simulation_time);
         NN               = length(d_T);
         ExcitationM      = zeros(NN,2);
         ExcitationM(:,1) = d_T;
         
         f                = f1+(f2-f1)/T.*d_T;  
         u_t              = 0.1 * sin( 2 *pi.*f.*d_T);
         ExcitationM(:,2) = u_t;
         simOut           = sim ('Simulator_Single_Axis1');

         % Input Controllato
         Mt = simOut.Mtot(1:end);

         % Output
         ax   = simOut.ax(1:end);
         q    = simOut.q(1:end);
         theta= simOut.theta(1:end);
         figure(1)
         plot(d_T,theta)
         grid on
         figure(2)
         plot(d_T,Mt)
         grid on
      Definizione sistema idgrey
parametres = {'Xu',Xu_MC;'Xq',Xq_MC;'Mu',Mu_MC;'Mq',Mq_MC;'Xd',Xd_MC;'Md',Md_MC;};

odefun = @quad2grey;
fcn_type = 'c';
sys0 = idgrey(odefun,parametres,fcn_type);   
         % Dataset Definition
       data = iddata( [q,ax],Mt,ts);
data.InputName = 'Input';
data.OutputName = {'q','ax'};

%trasformazione Dataset da time-domain to frequency domain
DATA_f = fft(data);

% Grey-box estimation
OPT = greyestOptions('SearchMethod','auto');
sys = greyest(DATA_f,sys0,OPT);

UNC_M  =  getcov(sys);
tr     =  trace(UNC_M);


