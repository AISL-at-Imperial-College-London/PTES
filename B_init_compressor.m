% Load the final state if it exists

state_file = 'x_comp.mat';
modelName = 'B_Comp';
try
    steady_state=x_comp;
catch
    if isfile(state_file)
        load(state_file, 'x_comp');  % Load the steady-state values
        % disp('Loaded steady_state for initialization.');
    else
        error('State file (x_comp.mat) not found! Run a steady-state simulation first.');
    end
end

load('splines.mat');
load('external.mat');

Param_COMP.arch = computer;
Param_COMP.splines= splines;

try
    Kp_tau;
    Ki_tau;
    stepping_time;
    DeltaW;
catch
    Kp_tau=1;
    Ki_tau=0.1;
    stepping_time=input('Step time:');
    DeltaW=input('Power change: ');
end