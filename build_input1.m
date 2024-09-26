function [covScore,INPUT,closed_loop_IO_data,data,sys,invM_] = build_input1(var_in,theta) 

        f_min = var_in(1);
        f_max = var_in(2);
        sim_time = var_in(3);
        simulation_time =sim_time;
        amp = 0.1;
        phi = 0;
        T = sim_time;
        ts = 0.004;
        t_nx = (0:ts:sim_time)';
        fun = @(t) amp*sin( 2*pi* ( ( f_max - f_min ) / sim_time.*t + f_min ).*t + phi );

        s = fun(t_nx);
        
 % reLOAD PARAMETERS(in model workspace)
% TO BE REPEATED AT EVERY STEP, AS GA MESSES THINGS UP!

model = 'Simulator_Single_Axis_mod1';
open_system(model);
MWS = get_param(model,'ModelWorkspace');
MWS.DataSource = 'MATLAB File';
MWS.FileName = 'DATA';
reload(MWS)

% Create structure INPUT to feed variables to simulink,
% after exporting to base workspace

% model = 'Simulator_Single_Axis_mod';
% open_system(model);
load_system(model);
INPUT.amp = amp;
INPUT.simulation_time=t_nx(end)-t_nx(1);
INPUT.ExcitationM = [t_nx,s];

A = [ theta(1) theta(2) 0;
      theta(3) theta(4) 0;
      0 1 0];
B = [theta(5);
    theta(6);
    0];
C = [1 0 0;
    0 1 0;
    0 0 1;
    theta(1) theta(2) 0];
D = [0 0 0 theta(5)]';

MWS = get_param(model,'ModelWorkspace');
assignin(MWS,'INPUT',INPUT)
assignin(MWS,'A',A)
assignin(MWS,'B',B)
assignin(MWS,'C',C)
assignin(MWS,'D',D)



closed_loop_IO_data = sim(model);


% Signal processing
                                                  
% Input Pitch Moment δlon (u) time domain

u_value = closed_loop_IO_data.Mtot;                %  Amplitude

% Output Pitch rate q (y1) time domain

q_value = closed_loop_IO_data.q;                   %  Amplitude 

% Output Acceleration a_x (y2) time domain

acc_value = closed_loop_IO_data.ax;                %  Amplitude


data_TD = iddata([q_value acc_value],u_value,0.004,'ExperimentName','Closed-loop');
data = fft(data_TD);

% initial conditions (static)
 
Xu0 = theta(1);
Xq0 = theta(2);
Mu0 = theta(3);
Mq0 = theta(4);
Xd0 = theta(5);
Md0 = theta(6);
 
parameters = {'Xu',Xu0;'Xq',Xq0;'Mu',Mu0;'Mq',Mq0;'Xd',Xd0;'Md',Md0};

sysid = idgrey(@quad2grey,parameters,'c');

% Grey-box estimation
sysid.InputName = {sprintf('%s_{lon}',char(948))};
sysid.OutputName = {'q','ax'};
% sysid.InputUnit = {'u'};
sysid.OutputUnit = {'rad/s','m/s^2'};
data.InputName = {sprintf('%s_{lon}',char(948))};
% data.InputUnit = {'u'};
data.OutputName = {'q','ax'};
data.OutputUnit = {'rad/s','m/s^2'};

opt = greyestOptions('InitialState','zero','DisturbanceModel','none', ...
    'EnforceStability',0,'EstimateCovariance',true,'Display','off');
opt.SearchOptions.MaxIterations = 20;
% opt.SearchOptions.MaxIterations = 7;

sys = greyest(data,sysid,opt);

invM_ = sys.CovarianceMatrix;           % built in covariance of parameters

covScore = trace(invM_);

toc

end

