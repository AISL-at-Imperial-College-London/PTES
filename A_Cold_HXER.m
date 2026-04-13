 function [sys,x0,str,tss]=A_Cold_HXER(t,x,u,flag,Param)

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

function [sys,x0,str,tss] = mdlInitialSizes(t,x,u, Param)

% This handles initialization of the function.
% Call simsize of a sizes structure.
N=Param.CHX.N;
sizes = simsizes;
sizes.NumContStates  = N;    % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = N;    % outputs of model
sizes.NumInputs      = 6;     % inputs of model
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

n=Param.CHX.N/2;
% Outputs
%%%%%%%%%%%% Enthalpy Formulation
h_CO2=x(1:n);

h_W = x(n+1:2*n);

P_CO2       = u(3);       % CO2 pressure kPa
% 
% T_W_out   =zeros(10,1);
% T_CO2_out  =zeros(10,1);

T_W_out = Param.water.T(h_W);
T_CO2_out = Param.splines.T(repelem(P_CO2,n)',h_CO2);

% for i=1:10
%     T_W_out(i)=refpropm('T','H',h_W(i),'P',P_W,'water');
%     T_CO2_out(i)=refpropm('T','H',h_CO2(i),'P',P_CO2,'CO2');
% end
% 
% 
sys = [T_CO2_out;T_W_out];

%%%%%%%%%%%%
% sys=x;

% ******************************************
% Derivatives
% ******************************************

function sys = mdlDerivatives(t,x,u,Param)

%% Pre-define parameters

UA = 5.7038e6*0.9;                 % Heat tran. Co. * Area W/K
n = Param.CHX.N/2;                  % number of elements
%%
% Inputs

T_W_in      = u(1);       % inlet water temperature K
T_CO2_in    = u(2);       % inlet CO2 temperature K
P_CO2       = u(3);       % CO2 pressure kPa
mCO2        = u(4);       % CO2 mass flow kg/s 
mW          = max(0,u(5));       % water mass flow kg/s
idle        = u(6);


water=Param.water;

%%%%%%%%%%%%%%% Enthalpy Formulation
rho_spline = Param.splines.rho;
h_spline = Param.splines.H;
T_spline = Param.splines.T;

hCO2in   = h_spline(P_CO2,T_CO2_in);
hWin     = water.h(T_W_in);

% hCO2in      =refpropm('H','T',T_CO2_in,'P',P_CO2,'CO2');
% hWin        =refpropm('H','T',T_W_in,'P',P_W,'water');

h_CO2   = x(1:n);
h_W     = x(n+1:2*n);

T_W = water.T(h_W);
rho_W = water.rho(T_W);

T_CO2 = T_spline(repelem(P_CO2,n)',h_CO2);
rho_CO2 = rho_spline(repelem(P_CO2,n)',T_CO2);

calc_HTC = Param.calc_HTC;

if calc_HTC == 1
    HTC_water=HTC(T_W,0,mW,rho_W,Param);
    HTC_CO2 = HTC(T_CO2,repelem(P_CO2,n)',mCO2,rho_CO2,Param);
    %Water side
    
    U=1./(1./HTC_CO2+1./flip(HTC_water));
else
    U = repelem(3.434397442853864e+06/n,n);
end

% disp('CHX Utot = ')
U_tot = sum(U);
% disp(U_tot)
% disp(t)

sys=zeros(n,1);
A=9.616;
V=A/2500;
if t >= 3450 && calc_HTC == 1
    CHX.U = U;
    CHX.Ut = U_tot;
    CHX.V = V;
    CHX.water=HTC_water;
    CHX.CO2 = HTC_CO2;
    
    assignin('base','CHX_Dat',CHX)
end
for i=1:n
    
    h_i=h_W(i);

    if i==1
        h_in=hWin;
    else
        h_in=h_W(i-1);
    end

    % T_i = refpropm('T','H',h_i,'P',P_W,'water');
    % rho_i=refpropm('D','H',h_i,'P',P_W,'water');

    T_i = T_W(i);
    rho_i=rho_W(i);

    j_cold=n+1-i;
    T_j = T_CO2(j_cold);
   
    segVol=V/2/n;

    mSeg=rho_i*segVol;

    dHdt=mW*(h_in-h_i)...
        -(U(j_cold)*A)/n*(T_i-T_j);

    sys(i+n)=dHdt/mSeg;
end

% CO2 Side

for j=1:n
    h_j=h_CO2(j);

    if j==1
        h_in=hCO2in;
    else
        h_in=h_CO2(j-1);
    end

    % T_j = refpropm('T','H',h_j,'P',P_CO2,'CO2');
    % rho_j = refpropm('D','H',h_j,'P',P_CO2,'CO2');

    T_j = T_CO2(j);
    rho_j=rho_CO2(j);

    i_hot = n+1-j;
    T_i = T_W(i_hot);

    % h_i = h_W(i_hot);
    % T_i = refpropm('T','H',h_i,'P',P_W,'water');



    segVolC=V/2/n;
    mSegC=rho_j*segVolC;

    dHdt=mCO2*(h_in-h_j)...
        +(U(j)*A/n)*(T_i-T_j);

    sys(j)=dHdt/mSegC;
end

% if t >= 633
%     disp("Time:")
%     disp(t)
%     disp("mCO2: ")
%     disp(mCO2)
%     disp("mW: ")
%     disp(mW)
%     disp(sys')
% end

function h=HTC(T,P,m,rho,params)
%%%%%%%%%%%%%%%%
% Calculates the heat transfer coefficient of a given flow
% Assumes semi-circular channels
% Based on Son et. al. 2020 doi: https://doi.org/10.1115/1.4047888
%%%%%%%%%%%%%%%%
if any(~isfinite(T)) || any(~isfinite(P)) || any(~isfinite(m)) || any(~isfinite(rho))
    error('HTC: non-finite inputs');
end
N=1000;
isWater = isscalar(P) && P==0;
if isWater
    v=params.water.nu(T);
    k=params.water.k(T);
    a=k./rho./params.water.Cp(T);
    R=1/1000;
else
    v=params.splines.nu(P,T);
    k=params.splines.k(P,T);
    a=k./rho./params.splines.Cp(P,T);

    R=0.5/1000;
end

L=0.67.*R;
Aflow=0.5*pi()*R^2;

Pr=v./a;

V_t=m./rho;
V_ch=V_t./N;

u=V_ch./Aflow;

Re=u.*L./v;

Nu=0.8405.*Re.^0.5704.*Pr.^1.08;

h=Nu.*k./L;

if any(~isfinite([v;k;a])) || any(v<=0) || any(rho<=0)
    error('HTC: bad property values (v, k, a, rho)');
end