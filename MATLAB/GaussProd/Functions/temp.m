function y = temp(x,coefs,shifts)

alpha = 13/30;

f = @(x,c,s) c * exp(-alpha*(x-s).^2);

y = zeros(size(x));

for k = 1:length(coefs)

    y = y + f(x,coefs(k),shifts(k));
    
end


    
    