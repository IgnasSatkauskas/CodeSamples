%example 2111 : Gaussian m2s1 times Gaussian m1s1
clear all; clc;
%close all;

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

format long

%alpha = 1/5;
alpha = 0.25;


% f(x) distribution can be arbitrary 
    
% f(x) is Gaussian (should have larger mean/sigma ratio)
mu_x = 2;
s_x = 1;
f = @(x) ( 1 / (s_x*sqrt(2*pi)) ) * exp( -(x - mu_x).^2 / (2*s_x^2) );


% g(y) distribution is Gaussian 
mu_y = 1;
s_y = 1;

savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha1over5/Examples/';

savename=['gauss_m2s1Xgauss_m1s1']

save_run = 1;

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
% savename

%% Moments 
moment0 = 1
moment1 = mu_x * mu_y
moment2 = mu_x^2*s_y^2 + mu_y^2*s_x^2 + s_x^2*s_y^2 


moments(W_kj, j_start, j_end, k_start, k_end, alpha)



%PLOTTING-----------------------------------------------------------------

%% FACTORS
% f(x) is Cauchy
t = linspace(-20,20,10000);
figure; plot(t,f(t),'k'); 
axis([-20 20 0 .5])

% g(x) is Gaussian 
g = @(x) (1/(sqrt(2*pi)*s_y)) * exp(-(x-mu_y).^2 / (2*s_y^2));
hold on; plot(t, g(t),'k'); 

%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-20 -10 0 10 20]; ax.YTick = [0 .1 .2 .3 .4 .5];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% PRODUCT
tt = linspace(-20,20,10000);

p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p2,'k');
axis([-20 20 0 .5])

% %setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-20 -10 0 10 20]; ax.YTick = [0 .1 .2 .3 .4 .5];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% SINGULARITY interval 2e-10 
tt = linspace(-2e-10,2e-10,5000);

p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);


figure; plot(tt,p2,'k');

axis([-2e-10 2e-10 0.75 1])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-2e-10 -1e-10 0 1e-10 2e-10]; ax.YTick = [.8 .9 1];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;


%% SINGULARITY interval 1e-26 
tt = linspace(-1e-26,1e-26,5000);

p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);


figure; plot(tt,p2,'k');

axis([-1e-26 1e-26 1.7 2])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-1 -.5 0 .5 1]*1e-26; ax.YTick = [1.7 1.8 1.9 2];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;


%% COEFFICIENTS
%old coefs
figure;
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj+1e-32),-32:.25:-2); %title('GMR matrix of coefficients (log scale)')
set(h,'lineStyle','none');

c=colorbar;
c.FontSize = 18;
c.Ticks=[-30:4:-2];

%colormap gray

%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])

%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-100 -50 0 50 100]; ax.YTick = [0 30 60 90];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
