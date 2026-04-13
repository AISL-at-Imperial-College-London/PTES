function T = extract_outlet_temps(u)

n=length(u)/2;

u1=u(1:n);
u2=u(n+1:2*n);

T1=u1(end);
T2=u2(end);
T=[T1,T2];
end