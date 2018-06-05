function y = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,t,alpha)

two = 2;

gscaling = @(x) sqrt(alpha/pi) * exp(-alpha*x^2);

np = length(t);
y = zeros(size(t));

epss = 1.d-15;
iwidth = ceil( sqrt(-log(epss)/alpha+log(two*alpha/pi)/alpha/4) );

for m=1:np
    y(m) = 0;
    for j=j_start:j_end
        fact = 2^(j/two);
        arg = t(m)*fact*fact;
        kk = ceil(arg);
        for k = max(k_start,kk-iwidth):min(k_end,kk+iwidth)
            y(m) = y(m) + W_kj(j-j_start+1, k-k_start+1) * gscaling(arg-k)*fact;
        end
    end
end



% !
% !
%     do m=1,np
%        gmval(m) = 0
%            do j=jsta,jsto
%               fact = two**(j/two)
%               arg  = xval(m)*fact*fact
%               kk = ceiling(arg)
% !
% ! takes the essential support of gscaling into account
% !
%               do k= max(ksta,kk-iwidth),min(ksto,kk+iwidth)
%               gmval(m) =gmval(m)+ cmatrix(k,j)*gscaling(arg-k,alpha)*fact
%               enddo
%            enddo
%      enddo       
