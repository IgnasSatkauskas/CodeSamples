function y = eval_GMR(W_vec,t,alpha)
%evaluate GMR representation given as vector of 
%significant coeficients

%Gaussian multiresolution basis element
% phi = @(k,j,t) 2^(j/2) * sqrt(alpha/pi) * exp( -alpha*(2^j*t - k).^2);
% j comes in W(*,1)
% k comes in W(*,2)
% weight is  W(*,3)

%alpha = 0.2;
const = sqrt(alpha/pi);

y = zeros(size(t));
for n = 1:length(W_vec(:,1))
    
     %y = y + W_vec(n,3) * phi(W_vec(n,2),W_vec(n,1),t);
     
     y = y + W_vec(n,3) * const * 2^(W_vec(n,1)/2) * exp( -alpha * (2^W_vec(n,1) * t - W_vec(n,2)).^2 ); 
     
end

%y = y.';




% %% scratch
% 
% A = [1 0 3; 4 5 0; 7 8 9; 10 11 0]
% 
% size(A);
% 
% ind = find(A<0.2)
% 
% x = mod(ind,size(A,1)) + 1;
% y = rem(ind,size(A,1));
% 
% 
% %lin_ind = sub2ind(size(A),[1 2 4],[2 3 3])
% 
% 
% [I, J] = ind2sub(size(A),ind);
% 
% [I, J]
% 
