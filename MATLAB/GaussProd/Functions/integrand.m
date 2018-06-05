function y = integrand(tau,f,s_y,mu_y,alpha,j,k)

%c1 = 1 / (sqrt(2*alpha)*s_y) ; % constant repeats in integrand definition 4 times; precompute
%integrand = @(tau) ( f( c1*2^(2-j)*2.^(-tau) ) .* exp( ( (-alpha*4.^(tau-2)) ./ (1 - 4.^(tau-2)) ) .* (k - c1*mu_y * 2.^(-tau+2) ).^2 ) ...
%                         + f( -c1*2^(2-j)*2.^(-tau) )  .* exp( ( (-alpha*4.^(tau-2)) ./ (1 - 4.^(tau-2)) ) .* (k + c1*mu_y * 2.^(-tau+2) ).^2 ) ) ./ ( sqrt(1-4.^(tau-2)) );

%alpha = 0.2;
%precompute repeated values
c1 = 1 / (sqrt(2*alpha)*s_y);
x_f = c1*2^(2-j)*2.^(-tau);
x_exp = (-alpha*4.^(tau-2)) ./ (1 - 4.^(tau-2));
x_sq = c1*mu_y * 2.^(-tau+2);

y = ( f(x_f) .* exp( x_exp .* (k - x_sq).^2 ) + f(-x_f) .* exp( x_exp .* (k + x_sq).^2 ) ) ./ ( sqrt(1-4.^(tau-2)) );



% % testing function handle speed
% 
% f = @(x) exp(x) .* sqrt(x);
% t = linspace(1,100,1000).';
% 
% 
% tic
% p = zeros(length(t),1);
% for i=1:100000
%     
%     p = p + f(t);
%     
% end
% toc
% 
% 
% tic
% p = zeros(length(t),1);
% for i=1:100000
%     
%     p = p + exp(t) .* sqrt(t);
%     
% end
%toc


