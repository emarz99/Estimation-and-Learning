%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANT-X SIMULATOR - MAIN                                                  %
% Authors:  Mattia Giurato (mattia.giurato@polimi.it)                     %
%           Paolo Gattazzo (paolo.gattazzo@polimi.it)                     %
% Date: 13/12/2017                                                        %
% Adapted to ANT-X 2DoF by:  Salvatore Meraglia (salvatore.meraglia@polimi.it)%
% Date: 22/12/2022                                                        %
%
% Further modified to include structure three-state identified longitudinal model
% 06/01/23 ML
% 
% Further modified to pick outputs with measurement error
% 03/01/24 ML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clearvars;
% close all;
% clear all;
addpath('datasets','common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');
% clc;
% Ts = 0.004
%% Model parameters

% Initial model (state: longitudinal velocity, pitch rate, pitch angle; input: normalised pitching moment; outputs: state and longitudinal acceleration)

Xu=-0.1068;

Xq=0.1192;

Mu=-5.9755;

Mq=-2.6478;

Xd=-10.1647;

Md=450.71;

A=[Xu, Xq, -9.81;
    Mu, Mq, 0;
    0, 1, 0];

B=[Xd;
    Md;
    0];

C=[1, 0, 0;
    0, 1, 0;
    0, 0, 1;
    Xu, Xq, 0]; 

D=[0;
    0
    0; 
    Xd];
%% Definizione sistema iniziale

% Noise

%noise.Enabler = 0;
noise.Enabler = 1;

noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]

noise.vel_stand_dev = noise.Enabler * 0.01;                               %[m/s]

noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);                   %[rad/s]

% Delays

delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

%% Load controller parameters

parameters_controller                    

%% M injection example (sweeep: first column time vector, secondo column time history of pitching moment) 

load ExcitationM

SetPoint=[0,0];


%% Launch SIMULATOR
n = 26646; %numero campioni
sample_time = 0.004;
simulation_time=106.5840;
 simOut = sim ('Simulator_Single_Axis1')



%% Mt : Momento di pitching totale
Mt = simOut.Mtot(2:end);

% Output
ax   = simOut.ax(2:end);
q    = simOut.q(2:end);
th   = simOut.theta1(2:end);
u    = simOut.u(2:end);



% Discretizzazione tempo
d_T = linspace(0,106.5840,n);
simulation_time = 106.5840;
ts = simulation_time/n;
fs = 1/ts;

plot(d_T,q)
% Definizione sistema idgrey


parametres = {'Xu',Xu;'Xq',Xq;'Mu',Mu;'Mq',Mq;'Xd',Xd;'Md',Md;};

odefun = @quad2grey;
fcn_type = 'c';
sys0 = idgrey(odefun,parametres,fcn_type);

%% Graphs
figure()
plot(d_T,Mt,'-k')
grid on
xlabel('Simulation Time [s]')
ylabel('Normalized Pitching Moment')
title('Normalized Pitching Moment')

figure()
plot(d_T,q,'-b')
grid on
xlabel('Simulation Time[s]')
ylabel('q [rad/s]')
title('Angular velocity q')

figure()
plot(d_T,ax,'-r')
grid on
xlabel('Simulation Time[s]')
ylabel('ax [m/s^2]')
title('Longitudinal acceleration')



%% Delete temporary files

if exist('slprj','dir')
    rmdir('slprj', 's')                                                    
end


%% Dataset definizione

data = iddata( [q,ax],Mt,ts);
data.InputName = 'Pitching Moment normalized';
data.OutputName = {'q','ax'};

%trasformazione Dataset da time-domain to frequency domain
DATA_f = fft(data);


%% Grey-box estimation
OPT = greyestOptions('SearchMethod','auto','EstimateCovariance',true);
sys = greyest(DATA_f,sys0,OPT)

opt = compareOptions('InitialCondition','zero');
compare(DATA_f,sys,Inf,opt);


%Computations of the matrix of covariances

UNC_M = getcov(sys);
T = sum(diag(UNC_M));
VETT = getpvec(sys);

% Computazione varianze in percentuale
for i = 1:6
a(i) = abs(100*UNC_M(i,i)/VETT(i));
end



%% Uncertainty analysis
UNC_Xu = sqrt(a(1))*3;
UNC_Xq = sqrt(a(2))*3;
UNC_Mu = sqrt(a(3))*3;
UNC_Mq = sqrt(a(4))*3;
UNC_Xd = sqrt(a(5))*3;
UNC_Md = sqrt(a(6))*3;

Xu_u = ureal('Xu',Xu,'Perc',UNC_Xu);
Xq_u = ureal('Xq',Xq,'Perc',UNC_Xq);
Mu_u = ureal('Mu',Mu,'Perc',UNC_Mu);
Mq_u = ureal('Mq',Mq,'Perc',UNC_Mq);
Xd_u = ureal('Xd',Xd,'Perc',UNC_Xd);
Md_u = ureal('Md',Md,'Perc',UNC_Md);
%%
A_u = [Xu_u, Xq_u, -9.81;
    Mu_u, Mq_u, 0;
    0, 1, 0];

B_u = [Xd_u;
      Md_u;
      0];

C_u = [1, 0, 0;
    0, 1, 0;
    0, 0, 1;
    Xu_u, Xq_u, 0]; 

D_u = [0;
       0;
       0; 
       Xd_u];
% Definizione sistema iniziale
IN_u = ss(A_u,B_u,C_u,D_u);

% Bode Diagram
figure(1)
bode(IN_u)
grid on
legend('Uncertain model','location','best')

figure(2)
pzplot(IN_u)
grid on
legend('Uncertain model','location','best')


