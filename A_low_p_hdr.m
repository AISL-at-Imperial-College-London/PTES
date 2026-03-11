function [sys,x0,str,tss]=A_low_p_hdr(t,x,u,flag,Param)

switch flag,

case 0,	% Initialize the states and sizes

   [sys,x0,str,tss] = mdlInitialSizes(t,x,u,Param);
     
case 3,   % Calculate the outputs
   
   sys = mdlOutputs(t,x,u,Param);

case 1,	% Obtain derivatives of states
   
   sys = mdlDerivatives(t,x,u,Param);

otherwise,

   sys = [];

end

% ******************************************
% Sub-routines or Functions
% ******************************************

% ******************************************
% Initialization
% ******************************************
function [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);

% This handles initialization of the function.
% Call simsize of a sizes structure.
sizes = simsizes;
sizes.NumContStates  = 1;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 1;     % outputs
sizes.NumInputs      = 3;     % inputs 
sizes.DirFeedthrough = 0;     % System is non-causal
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = 0;                   % Initialize the continuous states.
str = [];                     % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].

% ******************************************************************************************************************************
%  Outputs
% ******************************************************************************************************************************

function sys = mdlOutputs(t,x,u,Param);

% Outputs:

% States

sys(1) = x(1);

% ******************************************************************************************************************************
% Derivatives
% ******************************************************************************************************************************

function sys = mdlDerivatives(t,x,u,Param)

R = 188.9;            % gas constant for s CO2 (J/(kg·K))
%T = 420 + 273.15;     % High pressure side, use the compressor outlet T (K)
V = 80*3;               % Volume m3

% States

ph  = x(1);

% Inputs

dqRB = u(1); % m_comp kg/s
dqPB = u(2); % m_out from the compressor kg/s
T    = u(3); % temperature K

% Derivatives:
% sys(1)is dph/dt 
sys(1) = R*T/V*(dqRB-dqPB);      % In_pressure of the compressor Pa
%sys(1) = 5.581157407407407e+004/3*1e-6*(dqRB-dqPB);      % Out_pressure of the compressor Pa