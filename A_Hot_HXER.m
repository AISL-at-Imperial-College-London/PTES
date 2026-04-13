function [sys,x0,str,tss]=A_Hot_HXER(t,x,u,flag,Param)

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
n=Param.HHX.N;
sizes = simsizes;
sizes.NumContStates  = n;    % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = n+1;     % outputs of model
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

% Outputs
n=max(size(x));
% sys=repelem(x,2);
% idle = u(6);
% if t == 0
%     idle =0;
% end
% assignin('base','T_HHX_CO2',x(1:n/2))
% assignin('base','T_HHX_salt',x(n/2+1:n))
% if idle
%     sys(1:n/2) = u(2);
%     sys(n/2+1:n) = u(1);
% else
%     sys(1:n)=x;
% end
sys(1:n) = x;
sys(n+1) = u(5);
% ******************************************
% Derivatives
% ******************************************

function sys = mdlDerivatives(t,x,u,Param)

%% Pre-define parameters

% 40% KNO3 1.865 g/cm^3 = 1865 kg/m^3
% 60% NaNO3 1.9 g/cm^3 = 1900 kg/m^3 [https://wenku.baidu.com/view/fcdb3c8bcc22bcd126ff0c30.html?_wkts_=1709764789132&needWelcomeRecommand=1]
% Mixture mole mass 91.4 g/mol = 0.0914 kg/mol


UA = 1.9979e6;                  % heat transfer coefficienct W/K
num = Param.HHX.N/2;                  % number of elements
% rho_salt = 1712.33;            % density kg/m^3
% Cp_salt  = 1.556e+03;        % heat capacity J/(kg.K) （134 J/(mol.K)）  [https://elib.dlr.de/143749/1/PropertyAnalysis_SQM-DLR_final%20Report_v2.1.pdf;http://saimm.org.za/Conferences/Slags2004/029_Kawakami.pdf]

Cp_spline = Param.splines.Cp;
rho_spline = Param.splines.rho;
%% 
% Inputs
T_salt_in       = u(1);       % inlet salt temperature K
T_CO2_in        = u(2);       % inlet CO2 temperature K
P_CO2           = u(3);       % inlet CO2 pressure kPa
mCO2            = u(4);       % CO2 mass flow kg/s 
mSalt           = max(0,u(5));
idle            = u(6);
% salt mass flow kg/s
calc_HTC = Param.calc_HTC;
% rho_i = 140;
T_CO2=x(1:num);
T_salt=x(num+1:2*num);

Cp=Cp_spline(repelem(P_CO2,num)',T_CO2);
% Cp_0=Cp_spline(P_CO2,T_CO2_in);

rho=rho_spline(repelem(P_CO2,num)',T_CO2);

Cp_salt=(1.5405+3.0924e-4.*(T_salt-273.15)).*1000;
rho_salt=(2.106-6.6795e-4.*(T_salt-273.15)).*1000;
if calc_HTC == 1   
    HTC_CO2 = HTC(T_CO2,repelem(P_CO2,num)',mCO2,rho,Param);
    HTC_salt = HTC(T_salt,0,mSalt,rho_salt,Param);
    
    U = 1./(1./HTC_CO2+1./flip(HTC_salt));
else
    U = repelem(5.766498277894857e+06/num,num);
end

% disp('HHX Utot = ')
U_tot= sum(U);
% disp(U_tot)



A = 8.4328;
V=A/2500;
sys=zeros(1,2*num);
% Hot Side
if t>3450 && calc_HTC == 1
    HHX.U = U;
    HHX.Ut = U_tot;
    HHX.V = V;
    HHX.salt=HTC_salt;
    HHX.CO2 = HTC_CO2;
    assignin('base','HHX_Dat',HHX)
end

for i=1:num
    
    T_i=T_CO2(i);

    if i==1
        T_in=T_CO2_in;
    else
        T_in=T_CO2(i-1);
    end
    

    j_cold=num+1-i;

    T_j=T_salt(j_cold);
    
    sys(i)=(mCO2*Cp(i)*(T_in-T_i)...
        -(U(i)*A/num)*(T_i-T_j))/(rho(i)*V/2/num*Cp(i));
    % sys(i)=(mCO2*Cp(i)*(T_in-T_i)...
    %     -(UA/num)*(T_i-T_j))/(rho(i)*V/2/num*Cp(i));

end
    
for j=1:num
    
    T_j=T_salt(j);
    
    if j==1
        T_in=T_salt_in;
    else
        T_in=T_salt(j-1);
    end

    i_hot=num+1-j;
    
    T_i=T_CO2(i_hot);

    sys(num+j)=(mSalt*Cp_salt(j)*(T_in-T_j)...
                +(U(i_hot)*A/num)*(T_i-T_j))/(rho_salt(j)*V/2/num*Cp_salt(j));
    % sys(num+j)=(mSalt*Cp_salt(j)*(T_in-T_j)...
    %             +(UA/num)*(T_i-T_j))/(rho_salt(j)*V/2/num*Cp_salt(j));
    

end


function h=HTC(T,P,m,rho,params)
%%%%%%%%%%%%%%%%
% Calculates the heat transfer coefficient of a given flow
% Assumes semi-circular channels
% Based on Son et. al. 2020 doi: https://doi.org/10.1115/1.4047888
%%%%%%%%%%%%%%%%

N=1000;
if P==0
    mu=5.4047619e-8.*T.^2-8.0467958e-5.*T+3.1462e-2;
    v=mu./rho;
    k=0.55;
    Cp_salt=(1.5405+3.0924e-4.*(T-273.15)).*1000;
    a=k./rho./Cp_salt;
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

u_m=V_ch./Aflow;

Re=u_m.*L./v;

Nu=0.8405.*Re.^0.5704.*Pr.^1.08;

h=Nu.*k./L;
