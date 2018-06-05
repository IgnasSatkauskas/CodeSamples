function y = fplus(x);

n = length(x);
y = zeros(1,n);
for i = 1:n
    if x(i) >= 0
        y(i) = x(i);
    else
        y(i) = 0;
    end
    
end
