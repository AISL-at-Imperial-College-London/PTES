function [sys,x0,str,tss]=A_high_p_hdr(t,x,u,flag,Param)

switch flag

case 0	% Initialize the states and sizes

   [sys,x0,str,tss] = mdlInitialSizes(t,x,u,Param);
     
case 3   % Calculate the outputs
   
   sys = mdlOutputs(t,x,u,Param);

case 1	% Obtain derivatives of states
   
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
sizes.NumContStates  = 1;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 1;     % outputs
sizes.NumInputs      = 4;     % inputs 
sizes.DirFeedthrough = 0;     % System is non-causal
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = 0;                   % Initialize the continuous states.
str = [];                     % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].

% ******************************************************************************************************************************
%  Outputs
% ******************************************************************************************************************************

function sys = mdlOutputs(t,x,u,Param)

sys(1) = x(1);

function sys = mdlDerivatives(t,x,u,Param)

R = 188.9;            % gas constant for s CO2 (J/(kg·K))
%T = 420 + 273.15;     % High pressure side, use the compressor outlet T (K)
V = 0.15;               % Volume m3

% States

ph  = x(1);

% Inputs

m_out = u(1); % m_comp kg/s
m_cycle = u(2); % m_out from the compressor kg/s
T    = u(3); % temperature (K)
m_vent = u(4);
% Derivatives:

sys(1) = R*T/V*(m_out-m_cycle-m_vent);      % In_pressure of the compressor Pa
% sys(1) = R*T/V*(dqRB-dqPB);      % In_pressure of the compressor Pa