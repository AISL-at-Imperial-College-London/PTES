function Qlty = z_Qlty_CO2(u1,u2)
disp(u1)
if contains(computer,'MACA64')
    Q=PropsSI('Phase','T',u1,'P',u2*1000,'CO2');
    if Q ~= 1 && Q ~= 3
        Qlty=999;
    else
        Qlty=1;
    end
else
    Qlty = refpropm('Q','T',u1,'P',u2,'CO2');
end