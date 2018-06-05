function y = Uniform(a,b,x)
%uniform distribution on (a,b)


for i = 1:length(x)
    if (x(i) < b) && (x(i)>a)
        x(i) = 1/(b-a);
    else
        x(i) = 0;
    end
end
y=x;