%example Laplace times Gaussian
clear all; clc;
%close all;

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

format long

%alpha =1/5;
alpha = 0.25;

% f(x) distribution can be arbitrary 
    
% f(x) is a Laplace distribution
mu = 2;
b = 1;
f = @(x) (1/(2*b)) * exp(-abs(x-mu)/b);
% Laplace mean and std. deviation
mu_x = mu;
s_x = sqrt(2)*b;

% g(y) distribution is Gaussian 
mu_y = 1;
s_y = 1;

savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha1over5/Examples/';

%filename for saving runs

savename=['gauss_mneg2s1Xlaplace_m2b1']

save_run = 0;

% compute w_kj
k_start = -100;
k_end = 100;
j_start = -10;
j_end = 100;

W_kj = zeros(j_end-j_start+1,k_end-k_start+1);



tol_integral =10^-14;
%tol_weight = 10^-20;

tic
c1 = 1 / (sqrt(2*alpha)*s_y) ; % constant repeats in integrand definition 4 times; precompute
j_i = 0; %W_kj row index

% total number of funct evaluations
S = 0;
    
    
for j = j_start:j_end
    j
    j_i = j_i+1;
    k_i = 0; % W_kj column index
    for k = k_start:k_end
        k_i = k_i+1;
        %integrand(tau,f,s_y,mu_y,j,k)
        integ = @(tau) integrand(tau,f,s_y,mu_y,alpha,j,k);
        %AdaptiveIntegrator
        w_kj = ( ( 2^(-j/2)*log(2) ) / (sqrt(2*pi)*s_y) ) * AdaptiveIntegrator(integ,0,1,tol_integral);
        
        %count total function evaluations 
        [y,s] = AdaptiveIntegrator(integ,0,1,tol_integral);
        S = S + s;
     
%         if abs(w_kj) * 2^(j/2) < eps^2
%             w_kj=0;
%         end
        
        % zero out condition for pointwise error
        if j>0
            
            if abs(w_kj) * 2^(j/2) < eps*1e-3
                w_kj = 0;
            end
        else
            if abs(w_kj) * 2^(-j/2) < eps*1e-3
                w_kj =0;
            end
        end
        
        % put coefs into matrix
        W_kj(j_i,k_i) = w_kj;
        
    end
end


toc



% Make vector representation (only "significant" coefficients) 
tol_coef = 10^-16;
lin_ind = find(W_kj > tol_coef);
[j_i, k_i] = ind2sub(size(W_kj),lin_ind);

W_vec = zeros(length(lin_ind),3);
W_vec(:,1) = j_start + (j_i-1);
W_vec(:,2) = k_start + (k_i -1);
W_vec(:,3) = W_kj(lin_ind);


if save_run
    %save([savepath savename], 'W_vec')
    save([savepath savename], 'W_vec', 'W_kj', 'k_start','k_end','j_start','j_end')
    savename
end


% total integral check
[j_i_end,k_i_end] = size(W_kj);
total_integral =0;
for j_i = 1:j_i_end
    j = j_start + (j_i -1);
    total_integral = total_integral + sum(W_kj(j_i,:))*2^(-j/2);
end
total_integral
%savename

% total number of funct evaluations (one subdivision needs 20 evals)
number_func_eval = S*20

%% Moments 
moment0 = 1
moment1 = mu_x * mu_y
moment2 = mu_x^2*s_y^2 + mu_y^2*s_x^2 + s_x^2*s_y^2 


moments(W_kj, j_start, j_end, k_start, k_end, alpha)




%PLOTTING-----------------------------------------------------------------

%% FACTORS
% f(x) is Laplace
t = linspace(-20,20,10000);
figure; plot(t,f(t),'k'); 
axis([-20 20 0 .6])

% g(x) is Gaussian 
g = @(x) (1/(sqrt(2*pi)*s_y)) * exp(-(x-mu_y).^2 / (2*s_y^2));
hold on; plot(t, g(t),'k'); 

%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-20 -10 0 10 20]; ax.YTick = [0 .2 .4 .6];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% PRODUCT
tt = linspace(-20,20,10000);

%p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p,'k');

axis([-20 20 0 .5])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-20 -10 0 10 20]; ax.YTick = [.1 .3 .5];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% SINGULARITY
tt = linspace(-2e-10,2e-10,10000);

p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p,'k');

axis([-2e-10 2e-10 0.9 1.2])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-2e-10 -1e-10 0 1e-10 2e-10]; ax.YTick = [.9 1.0 1.1 1.2];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% SINGULARITY on 1e-28 d
tt = linspace(-1e-28,1e-28,5000);

%p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p2,'k');

axis([-1e-28 1e-28 2.24 2.36])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-1 -.5 0 .5 1]*1e-28; ax.YTick = [2.24 2.28 2.32 2.36];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% COEFFICIENTS

% %% add zero rows to the bottom to match Cauchy coefs i.e. j_start = -40
% j_start = -40;
% W_kj_add = zeros(j_end-j_start+1, k_end-k_start+1);
% W_kj_add(31:end,:) = W_kj; 
% W_kj = W_kj_add;
%%

figure;
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj+1e-32),-32:.25:-2); %title('GMR matrix of coefficients (log scale)')
set(h,'lineStyle','none');

c=colorbar;
c.FontSize = 18;
c.Ticks=[-30:4:-2];

colormap gray


%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])

%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-100 -50 0 50 100]; ax.YTick = [-30 0 30 60 90];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
