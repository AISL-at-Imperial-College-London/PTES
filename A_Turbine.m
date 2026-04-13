function [sys,x0,str,tss]=A_Turbine(t,x,u,flag,Param)

switch flag

case 0	  % Initialize the states and sizes
   
    [sys,x0,str,tss] = mdlInitialSizes(t,x,u,Param);

case 3   % Calculate the outputs
   
   sys = mdlOutputs(t,x,u,Param);
   
 case 1	% Obtain derivatives of states
% 
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

sizes = simsizes;
sizes.NumContStates  = 4;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 4;     % outputs of model    
sizes.NumInputs      = 6;     % inputs of model     
sizes.DirFeedthrough = 1;     % System is causal
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = 0;                   % Initialize the discrete states 
str = [];	                  % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].

% ******************************************
%  Outputs
% ******************************************
function sys = mdlOutputs(t,x,u,Param)


% Outputs

sys(1) = x(1);               % iso flowrate(kg/s)
sys(2) = x(2);        % iso outlet tempreature (K)
sys(3) = x(3);         % iso outlet quality (fraction)
sys(4) = x(4);    % Turbine work (W)


% ******************************************
% Derivatives
% ******************************************

function sys = mdlDerivatives(t,x,u,Param)

TauFlow = 0.01;
TauTemp = 0.01;     
TauQlty = 1   ;   
TauPow = 1;

% Inputs

T_inlet_iso    = u(1);      % iso inlet temperature (K)
P_inlet_iso    = u(2);      % iso inlet pressure (kPa)
Qlty_iso_inlet = u(3);      % iso inlet quality (fraction)
P_outlet_iso   = u(4);      % Turbine outlet pressure (kPa)
m_iso          = u(5);      % iso flowrate(kg/s)
eta            = u(6);

% Inlet quality

% if Qlty_iso_inlet < 1
%     error('Turbine inlet conditions are out of the operating range');
% else
Qlty_iso_inlet = 1;
% end

% Independent parameters

if contains(Param.arch,'MACA64')
    S_inlet_iso=PropsSI('S','T',T_inlet_iso,'P',P_inlet_iso*1000,'CO2');
    S_outlet_iso = S_inlet_iso;
    H_inlet_iso = PropsSI('H','T',T_inlet_iso,'P',P_inlet_iso*1000,'CO2');
    
    H_outlet_iso = PropsSI('H','P',P_outlet_iso*1000,'S',S_outlet_iso,'CO2');
    T_outlet_iso = PropsSI('T','P',P_outlet_iso*1000,'H',H_outlet_iso,'CO2');
    Turbine_work_isentropic_specific = H_inlet_iso - H_outlet_iso; 
    Turbine_work_actual_specific     = eta * Turbine_work_isentropic_specific;  % non-isentropic correction
    Cp_spline = Param.splines.Cp;
    Cp=Cp_spline(P_inlet_iso/1000,T_inlet_iso);
    Qlty_outlet                      = PropsSI('Phase','T',T_outlet_iso,'P',P_outlet_iso*1000,'CO2');

    % if Qlty_outlet == 1 || Qlty_outlet == 3
    %     error('Turbine outlet conditions are out of the operating range');
    % else
    %     Qlty_outlet = 1;
    % end
else
    S_inlet_iso   = refpropm('S','T',T_inlet_iso,'P',P_inlet_iso,'CO2');       % J/(kg.K)
    S_outlet_iso  = S_inlet_iso;                                               % isentropic expansion
    H_inlet_iso   = refpropm('H','T',T_inlet_iso,'S',S_inlet_iso,'CO2');       % J/kg

    % Isentropic expansion
    
    H_outlet_iso                     = refpropm('H','P',P_outlet_iso,'S',S_outlet_iso,'CO2');  % J/kg
    Turbine_work_actual_specific = H_inlet_iso - H_outlet_iso/eta;                             % J/kg
    % Turbine_work_actual_specific     = eta * Turbine_work_isentropic_specific;  % non-isentropic correction
    T_outlet_iso                     = refpropm('T','P',P_outlet_iso,'H',H_outlet_iso,'CO2');  % K
    Cp_spline = Param.splines.Cp;
    Cp=Cp_spline(P_inlet_iso/1000,T_inlet_iso);
    Qlty_outlet                      = refpropm('Q','T',T_outlet_iso,'P',P_outlet_iso,'CO2');  

    % if Qlty_outlet < 1
    %     error('Turbine outlet conditions are out of the operating range');
    % else
    %     Qlty_outlet = 1;
    % end
end


% % Isentropic expansion
% 
% S_outlet_gas_sat    = refpropm('S','P',P_outlet_iso,'Q',1,'CO2');             % gas saturation entropy J/(kg.K)
% S_outlet_liquid_sat = refpropm('S','P',P_outlet_iso,'Q',0,'CO2');             % liquid saturation entropy J/(kg.K)
% H_outlet_gas_sat    = refpropm('H','P',P_outlet_iso,'Q',1,'CO2');             % gas saturation enthalpy J/kg
% H_outlet_liquid_sat = refpropm('H','P',P_outlet_iso,'Q',0,'CO2');             % liquid saturation enthalpy J/kg
% 
% if S_outlet_iso < S_outlet_gas_sat
% 
%     Qlty_outlet                      = (S_outlet_iso - S_outlet_liquid_sat)/(S_outlet_gas_sat - S_outlet_liquid_sat);
%     H_outlet_iso                     = Qlty_outlet * H_outlet_gas_sat + (1-Qlty_outlet) * H_outlet_liquid_sat;        % J/kg
%     Turbine_work_isentropic_specific = H_inlet_iso - H_outlet_iso;                                                    % J/kg
%     Turbine_work_actual_specific     = eta * Turbine_work_isentropic_specific;                         % non-isentropic correction
% 
%     H_outlet_iso                     = H_inlet_iso - Turbine_work_actual_specific;                                    % J/kg
% 
%     if H_outlet_iso < H_outlet_gas_sat
%         Qlty_outlet  = (H_outlet_iso - H_outlet_liquid_sat)/(H_outlet_gas_sat - H_outlet_liquid_sat);
%         %T_outlet_iso = refpropm('T','P',P_outlet_iso,'X',Qlty_outlet,'CO2');
%         T_outlet_iso = refpropm('T','P',P_outlet_iso,'Q',0,'CO2');                                                    % ??????? is it correct that Q is 0? this can be a two-phase area
%     else
%         Qlty_outlet  = 1;
%         T_outlet_iso = refpropm('T','P',P_outlet_iso,'H',H_outlet_iso,'CO2');
%     end    
% else
%     Qlty_outlet = 1;
%     H_outlet_iso                     = refpropm('H','P',P_outlet_iso,'S',S_inlet_iso,'CO2');
%     Turbine_work_isentropic_specific = H_inlet_iso - H_outlet_iso;
%     Turbine_work_actual_specific     = eta * Turbine_work_isentropic_specific; % non-isentropic correction
%     H_outlet_iso                     = H_inlet_iso - Turbine_work_actual_specific;
%     T_outlet_iso                     = refpropm('T','P',P_outlet_iso,'H',H_outlet_iso,'CO2');
% end

kappa = 1.24;

iso_P_ratio = P_outlet_iso/P_inlet_iso;
Power = ((m_iso * Cp * T_inlet_iso)*eta * (1/iso_P_ratio^((kappa-1)/kappa)-1))/1000;

if m_iso<10
    Power_specific = 0;
else
    Power_specific = Power/m_iso*1000;
end
if contains(Param.arch,'MACA64')
    T_outlet_iso = PropsSI('T','P',P_outlet_iso*1000,'H',H_inlet_iso-Power_specific,'CO2');
else  
    T_outlet_iso = refpropm('T','P',P_outlet_iso,'H',H_inlet_iso-Power_specific,'CO2');  % K
end
Turbine_work_actual = Turbine_work_actual_specific * m_iso;

sys(1) = TauFlow*(m_iso - x(1)) ;              % iso flowrate(kg/s)
sys(2) = TauTemp*(T_outlet_iso - x(2)) ;       % iso outlet tempreature (K)
sys(3) = TauQlty*(Qlty_outlet - x(3)) ;        % iso outlet quality (fraction)
sys(4) = TauPow*(Power - x(4)) ;% Turbine work (W)

