function x=convertT_H(state,P,F)
n=max(size(state));
x=zeros(1,n);

if mean(state)<=1000
    for i=1:n
        x(i)=refpropm('H','T',state(i),'P',P,F);
    end
else
    for i=1:n
        x(i)=refpropm('T','H',state(i),'P',P,F);
    end

end
end
