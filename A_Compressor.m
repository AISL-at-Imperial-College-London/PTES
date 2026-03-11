function [sys,x0,str,tss]=A_Compressor(t,x,u,flag,Param)

%% Structure
switch flag

case 0	  % Initialize the states and sizes

   [sys,x0,str,tss] = mdlInitialSizes(t,x,u,Param);

case 3   % Calculate the outputs
 
   sys = mdlOutputs(t,x,u,Param);
   
case 1	  % Obtain derivatives of states

   sys = mdlDerivatives(t,x,u,Param);

    otherwise

   sys = [];

end

% ******************************************
% Sub-routines or Functions
% ******************************************
% ******************************************
% Initialization
% ******************************************

function [sys,x0,str,tss] = mdlInitialSizes(t,x,u,Param)

% This handles initialization of the function.
% Call simsize of a sizes structure.

sizes = simsizes;
sizes.NumContStates  = 6;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 12;    % outputs of model
sizes.NumInputs      = 7;     % inputs of model
sizes.DirFeedthrough = 1;     % System is causal
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = 0;                   % Initialize the states
str = [];	                  % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].

% ******************************************
%  Outputs
% ******************************************

function [sys] = mdlOutputs(t,x,u,Param)

global gpr_comp_map_PR gpr_comp_map_eta speed_mean speed_std flow_mean flow_std PR_mean PR_std eta_mean eta_std

% global m_in m_out m_rec_ss omega_scaled flow_scaled p_ratio_ss eta T_out Power surge_d p_comp torque_comp

Dt         = 2.2252;                       % impeller diameter m     （torque equations)  
sigma_slip = 0.9;                          % compressor slip

Valve_in_gain  = 0.1071*1.34;                                    % opening gain
Valve_out_gain = 0.0459*1.9;                                     % opening gain

% Inputs

torque_drive       = u(1); %
Inflow_opening     = u(2); %
Outflow_opening    = u(3); %
Recycle_opening    = u(4); %
In_pres            = u(5); %
Out_pres           = u(6); %
T_in               = u(7); %
% Cp                 = u(8); 

% States

p1         = x(1);  %
p2         = x(2);  %
m_comp     = x(3);  %
p_ratio    = x(4);  %
omega_comp = x(5);  %
m_rec      = x(6);  %

%% Algebraic equations

% Valves

Valve_rec_gain = Valve_out_gain/2;                                                     % opening gain
% m_in     = max(0,Valve_in_gain * Inflow_opening * sqrt(abs(In_pres - p1)) * (In_pres - p1)/abs(In_pres - p1));      % Inflow valve
% m_out    = max(0,Valve_out_gain * Outflow_opening * sqrt(abs(p2 - Out_pres)) * (p2 - Out_pres)/abs(p2 - Out_pres)); % Outflow valve
% % m_rec_ss = Valve_rec_gain * Recycle_opening * sqrt(abs(p2 - p1)) * (p2 - p1)/abs(p2 - p1);                   % Recycle valve

if omega_comp >= 0.1
    m_in     = max(0,Valve_in_gain * Inflow_opening * sqrt(abs(In_pres - p1)) * (In_pres - p1)/abs(In_pres - p1));      % Inflow valve
    m_out    = max(0,Valve_out_gain * Outflow_opening * sqrt(abs(p2 - Out_pres)) * (p2 - Out_pres)/abs(p2 - Out_pres)); % Outflow valve
    m_rec_ss = Valve_rec_gain * Recycle_opening * sqrt(abs(p2 - p1)) * (p2 - p1)/abs(p2 - p1);                   % Recycle valve
else
 m_in     = max(0,Valve_in_gain * Inflow_opening * sqrt(abs(In_pres - p1)) * (In_pres - p1)/abs(In_pres - p1));      % Inflow valve
    m_out    = max(0,Valve_out_gain * Outflow_opening * sqrt(abs(p2 - Out_pres)) * (p2 - Out_pres)/abs(p2 - Out_pres)); % Outflow valve
    m_rec_ss = Valve_rec_gain * Recycle_opening * sqrt(abs(p2 - p1)) * (p2 - p1)/abs(p2 - p1);                   % R
end
% compressor pressure ratio

omega_scaled = (omega_comp/2/pi*60 - speed_mean)/speed_std;                                 % omega -- rad/s to rpm
flow_scaled = (m_comp - flow_mean)/flow_std;
% p_ratio_ss = predict(gpr_comp_map_PR,[omega_scaled flow_scaled]) * PR_std + PR_mean *1.3;   % set an offset for p_ratio
% p_ratio_ss = max(p_ratio_ss, 1);

% ******* Isentropic efficiency ********

eta = predict(gpr_comp_map_eta,[omega_scaled flow_scaled]) * eta_std + eta_mean;
eta = max(eta, 0.4);
eta = min(eta, 0.9);

% ******* Outlet temperature  **********
kappa = 1.24;
T_out = T_in * (1 + ((p_ratio)^((kappa-1)/kappa) - 1)/eta);                     % Outlet temperature K

% **********      Power       **********
Cp_spline = Param.splines.Cp;
Cp=Cp_spline(p1/1000,T_in);

Power = ((m_comp * Cp * T_in)/eta * (p_ratio^((kappa-1)/kappa)-1)) /1000;       % Power (kW); m_comp (kg/s); Cp (J/(kg.K))

% compressor torque

torque_comp = sigma_slip*(Dt/2)^2*omega_comp*abs(m_comp); % (N.m) (rad/s) (kg/s)
%torque = 9548 * Power/(omega_comp/2/pi*60);

% Outputs

sys(1) = x(1)/1000;               % p1 kPa
sys(2) = x(2)/1000;               % p2 kPa
sys(3) = x(3);                    % m_comp kg/s
sys(4) = x(4);                    % p_ratio  
sys(5) = x(5);                    % omega_comp  rad/s
sys(6) = m_in;                    % m_in kg/s
sys(7) = m_out;                   % m_out kg/s
sys(8) = m_rec;                   % m_rec kg/s
sys(9) = Power;                   % Power kW
sys(10) = T_out;                  % T_out K
sys(11) = eta;                    % eta
sys(12) = torque_comp;            % torque_comp N.m

% ******************************************
% Derivatives
% ******************************************

function sys = mdlDerivatives(t,x,u,Param)

global gpr_comp_map_PR speed_mean speed_std flow_mean flow_std PR_mean PR_std 

% global p2 m_comp p_ratio m_rec 
% global m_in m_out m_rec_ss p_ratio_ss p_comp torque_comp

Dt              = 2.2252;                       % impeller diameter m     （torque equations)  
sigma_slip      = 0.9;                          % compressor slip

Valve_in_gain   = 0.1071*1.34;                                    % opening gain
Valve_out_gain  = 0.0459*1.9;                                     % opening gain

VolumeT1        = 20; % volume m^3  
VolumeT2        = 10; % volume m^3
AdivL           = 1.12147223819255e-05; % Duct cross section divided by length m

J               = 590.2839;                     % shaft inertia kg*m^2     (torque[N.m] = I[kg*m^2]*omega[rad/s])

tauComp    = 0.9;                          % time constant of the pressure ratio   
tauRecycle = 1;                            % time constant of the recycle        
SpeedSound = 404;                          % speed of sound m/s

torque_drive       = u(1); %
Inflow_opening     = u(2); %
Outflow_opening    = u(3); %
Recycle_opening    = u(4); %
In_pres            = u(5); %
Out_pres           = u(6); %




p1         = x(1);  %
p2         = x(2);  %
m_comp     = x(3);  %
p_ratio    = x(4);  %
omega_comp = x(5);  %
m_rec      = x(6);  %



Valve_rec_gain = Valve_out_gain/2;   % opening gain

m_in     = max(0,Valve_in_gain * Inflow_opening * sqrt(abs(In_pres - p1)) * (In_pres - p1)/abs(In_pres - p1));      % Inflow valve
m_out    = max(0,Valve_out_gain * Outflow_opening * sqrt(abs(p2 - Out_pres)) * (p2 - Out_pres)/abs(p2 - Out_pres)); % Outflow valve
m_rec_ss = Valve_rec_gain * Recycle_opening * sqrt(abs(p2 - p1)) * (p2 - p1)/abs(p2 - p1);                   % Recycle valve

% compressor pressure ratio

omega_scaled = (omega_comp/2/pi*60 - speed_mean)/speed_std;                                 % omega -- rad/s to rpm
flow_scaled = (m_comp - flow_mean)/flow_std;
p_ratio_ss = predict(gpr_comp_map_PR,[omega_scaled flow_scaled]) * PR_std + PR_mean *1.3;   % set an offset for p_ratio
% p_ratio_ss = max(p_ratio_ss, 1);

if omega_comp <205|| m_comp< 0.1
     p_ratio_ss = 0.6/205*(omega_comp)+1;
end
p_comp = p_ratio*p1;

% compressor torque

torque_comp = sigma_slip*(Dt/2)^2*omega_comp*(m_comp)+omega_comp*0.001;   % (N.m) (rad/s) (kg/s)



sys(1) = SpeedSound^2/VolumeT1*(m_in + m_rec - m_comp);     % Delta suction pressure
sys(2) = SpeedSound^2/VolumeT2*(m_comp - m_rec - m_out);    % Delta discharge pressure
sys(3) = 0.2*AdivL*(p_comp-p2);                                 % Delta mass in compressor sys(3) = 0.1*AdivL*(p_comp-p2)
sys(4) = tauComp*(p_ratio_ss-p_ratio);                      % Delta pressure ratio     
sys(5) = 1/J*(torque_drive-torque_comp);                    % Delta speed rad/s   
sys(6) = tauRecycle*(m_rec_ss - m_rec);                     % Delta recycle  


if omega_comp <0.5|| m_comp< 0.1
    sys(3) = -0.01*m_comp + torque_drive/10000;
    sys(5) = 1/J*(torque_drive-torque_comp)-x(5);  
end