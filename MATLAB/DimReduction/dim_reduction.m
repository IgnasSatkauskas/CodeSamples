function Rk = dim_reduction(Rn,k,option)

%input: Rn - matrix of vectors(as rows) in n dim
%       k - desired dim
%       oprion - (char) dim reduction option


[n,m] = size(Rn);
% n - number of points
% m = dimention of each point

%simple reduction
if option == 's' 
    
    ind = rand_pick(1,m,k);
    Rk = Rn(:,ind);

end

if option == 'p'
    
    %APPLY P

    q = min( 1, (log(n))^2 / m ); %prob of not zero entry in P (density)
    P = sprandn(k,m,q);%*sqrt(1/q);
    %full(P)
    %q
    
    Rk = ( P*Rn' )';
    
end



if option == 'phd'


    % APPLY D

    m2 = floor(m/2);
    ind = rand_pick(1,m,m2); %pick coordinates to multiply by -1
    Rn(:,ind) = Rn(:,ind)*(-1); %multiply selected coord by -1
    
    % APPLY H
    H=zeros(m,m);
    for i = 2:m+1
        for j=2:m+1
            ii = dec2bin(i);
            jj = dec2bin(j);
            dig = min(length(ii),length(jj));
            ii = ii(length(ii)-dig+1:end);
            jj = jj(length(jj)-dig+1:end);
            expon = ii*jj';
            H(i-1,j-1) = (1/sqrt(m)) * (-1)^expon;
        end
    end
    H 
   
    %APPLY P

    q = min( 1, (log(n))^2 / m ); %prob of not zero entry in P (density)
    P = sqrt(1/q)*sprandn(k,m,q);
    Rk = ( P*Rn' )';
    
end

