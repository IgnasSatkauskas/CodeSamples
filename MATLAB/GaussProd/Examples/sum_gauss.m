% for LaplaseXGumbel example:
%function g: laplace approximation as sum of Gaussians
% c_n coeficients of Gaussians: must be row
% s_n sigma's of Gaussians: must be col 
% x point(s) of evaluation 


function y = sum_gauss(c_n, s_n, mu, x)

y = zeros(size(x));

for n = 1:length(x)
    
    y(n) = c_n * exp(-(x(n)-mu)^2 ./ (2*s_n.^2) );
    
end

