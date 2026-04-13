function [h,Re]=nusselt(T,P,m,rho,params)
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

end