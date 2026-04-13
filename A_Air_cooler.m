function [sys,x0,str,tss]=A_Air_cooler(t,x,u,flag,Param)

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
N=Param.COL.N;
sizes = simsizes;
sizes.NumContStates  = N;    % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = N;     % outputs of model
sizes.NumInputs      = 5;     % inputs of model
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
n=Param.COL.N/2;
% Outputs
P_CO2           = u(3);       % CO2 pressure kPa

T_spline = Param.splines.T;

T_out = T_spline(repelem(P_CO2,n)',x(1:n));

T_air_out = Param.air.T(x(n+1:2*n));
% T_out=zeros(10,1);
% for i=1:10
%     T_out(i)=refpropm('T','H',x(i),'P',P_CO2,'CO2');
% end

sys = [T_out; T_air_out];

% sys=x;

% ******************************************
% Derivatives
% ******************************************

function sys = mdlDerivatives(t,x,u,Param)
%% Inputs

T_CO2_in        = u(1);       % inlet CO2 temperature K
T_air_in        = u(2);       % cold inlet CO2 temperature K
P_CO2           = u(3);       % CO2 pressure kPa
mass            = u(4);       % CO2 mass flow kg/s 
m_air           = u(5);

%% Params
h_spline = Param.splines.H;
T_spline = Param.splines.T;
rho_spline = Param.splines.rho;
Air = Param.air;

V = 10;                       % Volume m^3
UA = 500000*0.1111;                       % heat transfer coefficienct*Area W/K
n = Param.COL.N/2;                  % number of elements


%% States

h           = x(1:n);
h_air       = x(n+1:2*n);

hCO2in = h_spline(P_CO2,T_CO2_in);
hAin   = Air.h(T_air_in);


T = T_spline(repelem(P_CO2,n)',h);
rho = rho_spline(repelem(P_CO2,n)',T);

T_air = Air.T(h_air);
rho_air = Air.rho(T_air);

sys=zeros(n,1);
for i=1:n

    h_i=h(i);
    
    if i==1
        h_in=hCO2in;
    else
        h_in=h(i-1);
    end

    T_i = T(i);
    rho_i=rho(i);

    j_cold=n+1-i;
    T_j = T_air(j_cold);
    
    segVol=V/n;
    mSeg = rho_i*segVol;

    dHdt=mass*(h_in-h_i)-(UA/n)*(T_i-T_j);

    sys(i)=dHdt/mSeg;
end


for j=1:n
    h_j=h_air(j);

    if j==1
        h_in=hAin;
    else
        h_in=h_air(j-1);
    end
    
    T_j = T_air(j);
    rho_j = rho_air(j);

    i_hot=n+1-j;
    T_i=T(i_hot);

    segVolC=V/n;
    mSegC = rho_j*segVolC;

    dHdt=m_air*(h_in-h_j)...
        +(UA/n)*(T_i-T_j);

    sys(n+j)=dHdt/mSegC;
end