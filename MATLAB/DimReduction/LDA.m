function [M,s,S,Pi,drop_coor] = LDA(T, T_idx_g, del)

%del = .8;
[N_T,p] = size(T);
%Split T by classes

for i = 1:6
    idx = find(T_idx_g == i);
    eval(['T' int2str(i) '=T(idx,:);']); % Ti = T(idx,:);
end

%Find class centroids(M), and std's of each coordinate (for each class)  

M = zeros(6,p); % class centroids as rows of M
S = zeros(6,p);% ith row of S = std's of each coordinate for ith class
for i = 1:6 
    eval(['Tk = T' int2str(i) ';']); %Tk = Ti;
    N_k(i) = size(Tk,1); %size of each class
    M(i,:) = mean(Tk); % means of Tk columns (class centroids)
    S(i,:) = std(Tk); % std's of Tk columns
    
end

%pool within-class std's for each coordinate

S = S.^2;
s = (N_k - 1 ) * S;
s = s/(N_T-6); %pooled variances(std's squared) of each coordinate 1 by p

%Find class priors
Pi = N_k/N_T;

%Shrinking class centroids M

s0 = median(s);
m = mean(T); % mean of columns of T (overall means of each coord)
D = zeros(6,p);

for i = 1:6
    mk(i) = sqrt( 1/N_k(i) - 1/N_T );
    D(i,:) = (1/mk(i)) * ( M(i,:) - m ) ./ (s + s0);
    D(i,:) = sign( D(i,:) ) .* fplus( abs(D(i,:)) - del  );
    M(i,:) = m + mk(i) * (s-s0) .* D(i,:);  %shrunken centroids 
    
end


zero_D_col = find( sum( abs(D) ) == 0 ); % zero columns of D
drop_coor = length ( zero_D_col ); % number of coord droped




end

    
    
    