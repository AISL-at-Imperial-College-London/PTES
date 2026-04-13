function [sys,x0,str,tss]=A_REC_h(t,x,u,flag,Param)

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
n=Param.REC.N;
sizes = simsizes;
sizes.NumContStates  = n;    % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = n;    % outputs of model
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
n=Param.REC.N/2;
% Outputs
h_hot   =    x(1:n);
h_cold  =    x(n+1:2*n);

P_hot=  u(3);
P_cold= u(4);

%Convert enthalpy to temperature
T_spline=Param.splines.T;

T_hot_out=T_spline(repelem(P_hot,n)',h_hot);
T_cold_out=T_spline(repelem(P_cold,n)',h_cold);

sys = [T_hot_out;T_cold_out];

% ******************************************
% Derivatives
% ******************************************

function sys = mdlDerivatives(t,x,u,Param)

%% Pre-define parameters

% Old parameters
% V = 5.3 ;                        % Volume m^3
% UA = 9e6*1000;                      % heat transfer coefficienct*Area W/K

%% 
% Inputs

T_hot_in    = u(1);       % hot inlet salt temperature K
T_cold_in   = u(2);       % cold inlet CO2 temperature K
P_hot       = u(3);       % hot CO2 pressure kPa
P_cold      = u(4);       % cold CO2 pressure kPa
mass        = u(5);       % CO2 mass flow kg/s
idle        = u(6);

calc_HTC = Param.calc_HTC;

h_spline    = Param.splines.H;
rho_spline  = Param.splines.rho;
T_spline    = Param.splines.T;

hHot_in     = h_spline(P_hot,T_hot_in); 
hCold_in    = h_spline(P_cold,T_cold_in);

% hHot_in    = refpropm('H','T',T_hot_in,'P',P_hot,'CO2');
% hCold_in   = refpropm('H','T',T_cold_in,'P',P_cold,'CO2');

n=Param.REC.N/2;

h_hot       = x(1:n);
h_cold      = x(n+1:2*n);

T_hot       = T_spline(repelem(P_hot,n)',h_hot);
T_cold      = T_spline(repelem(P_cold,n)',h_cold);

% T_hot=zeros(1,10);
% T_cold=zeros(1,10);
% for i=1:10
%     T_hot(i)=refpropm('T','H',h_hot(i),'P',P_hot,'CO2');
%     T_cold(i)=refpropm('T','H',h_cold(i),'P',P_cold,'CO2');
% end
sys         = zeros(2*n,1);

rho_hot     = rho_spline(repelem(P_hot,n)',T_hot);
rho_cold    = rho_spline(repelem(P_cold,n)',T_cold);
if calc_HTC == 1
    HTC_hot = HTC(T_hot,repelem(P_hot,n)',mass,rho_hot,Param);
    HTC_cold= HTC(T_cold,repelem(P_cold,n)',mass,rho_cold,Param);
    
    U           = 1./(1./HTC_hot+1./flip(HTC_cold));
else
    U = repelem(4.405214457118595e6/n,n);
end
% disp('REC Utot = ')
U_tot       = sum(U);
% disp(U_tot)
A           = 11;
V           = A/2500; %specific volume for PCHEs 2500 m2/m3. Source Oh et. al. DOI: 10.1115/1.3126780
if t>3450 && calc_HTC == 1 
    REC.U = U;
    REC.Ut = U_tot;
    REC.V = V;
    REC.hot=HTC_hot;
    REC.cold = HTC_cold;
    
    assignin('base','REC_Dat',REC)
end
% Hot side
for i=1:n
    h_i = h_hot(i);

    if i==1
        h_in=hHot_in;
    else
        h_in=h_hot(i-1);
    end

    T_i = T_hot(i);
    rho_i=rho_hot(i);
    
    j_cold=n+1-i;
    T_j=T_cold(j_cold);

    segVolH=V/2/n;

    mSegH = rho_i*segVolH;

    dHdt=mass*(h_in-h_i)...
        -U(i)*A/n*(T_i-T_j);
    
    sys(i)=dHdt/mSegH;
end

% Cold side
for j=1:n
    h_j   = h_cold(j);
    
    if j==1
        h_in=hCold_in;
    else
        h_in=h_cold(j-1);
    end

    T_j = T_cold(j);
    rho_j = rho_cold(j);

    i_hot = n+1 - j;
    T_i = T_hot(i_hot);

    
    segVolC=V/2/n;

    mSegC=rho_j*segVolC;

    dHdt=mass*(h_in-h_j)...
        +(U(i_hot)*A/n)*(T_i-T_j);

    sys(n+j)=dHdt/mSegC;
end


function h=HTC(T,P,m,rho,params)
%%%%%%%%%%%%%%%%
% Calculates the heat transfer coefficient of a given flow
% Assumes semi-circular channels
% Based on Son et. al. 2020 doi: https://doi.org/10.1115/1.4047888
%%%%%%%%%%%%%%%%
R=0.5/1000;

L=0.67.*R;
Aflow=0.5*pi()*R^2;

N=1000;
% fluid=params.fluid_name;

% v=refpropm('$','T',T,'P',P,fluid)*1e-4;
% a=refpropm('%','T',T,'P',P,fluid)*1e-4;
% rho=refpropm('D','T',T,'P',P,fluid);
% k=refpropm('L','T',T,'P',P,fluid);
v=params.splines.nu(P,T);
% rho=params.splines.rho(P,T);
k=params.splines.k(P,T);
a=k./rho./params.splines.Cp(P,T);
Pr=v./a;

V_t=m./rho;
V_ch=V_t./N;

u=V_ch./Aflow;

Re=u.*L./v;

Nu=0.8405.*Re.^0.5704.*Pr.^1.08;

h=Nu.*k./L;