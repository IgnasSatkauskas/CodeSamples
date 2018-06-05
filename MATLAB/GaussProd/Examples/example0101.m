%example 0101 : Gaussian m0s1 times Gaussian m0s1
% Bessel K example
clear all; clc;
%close all;

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

format long

%alpha = 1/5;
%alpha = 0.2;
alpha = 0.25;


% f(x) distribution can be arbitrary 
    
% f(x) is Gaussian (should have larger mean over sigma ratio)
mu_x = 0;
s_x = 1;
f = @(x) ( 1 / (s_x*sqrt(2*pi)) ) * exp( -(x - mu_x).^2 / (2*s_x^2) );


% g(y) distribution is Gaussian 
mu_y = 0;
s_y = 1;

savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha1over5/Examples/';

savename=['gauss_m0s1Xgauss_m0s1']

save_run = 0;

% compute w_kj
k_start = -100;
k_end = 100;
j_start = -10;
j_end = 100;

W_kj = zeros(j_end-j_start+1,k_end-k_start+1);


tol_integral =10^-12;
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
            
            if abs(w_kj) * 2^(j/2) < eps*1e-6
                w_kj = 0;
            end
        else
            if abs(w_kj) * 2^(-j/2) < eps*1e-6
                w_kj =0;
            end
        end
        
        % put coefs into matrix
        W_kj(j_i,k_i) = w_kj;
        
    end
end


toc



% Make vector representation (only "significant" coefficients) 
%tol_coef = 1e-38;
lin_ind = find(W_kj ~= 0);
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
% f(x) is Gaussian
t = linspace(-10,10,10000);
figure; plot(t,f(t),'k'); 
axis([-10 10 0 .5])

% g(x) is Gaussian 
g = @(x) (1/(sqrt(2*pi)*s_y)) * exp(-(x-mu_y).^2 / (2*s_y^2));
hold on; plot(t, g(t),'k'); 

%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-10 -5 0 5 10]; ax.YTick = [0 .1 .2 .3 .4 .5];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% PRODUCT
tt = linspace(-10,10,5000);

p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);


figure; plot(tt,p2,'k');
axis([-10 10 0 2.5])

% %setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-10 -5 0 5 10]; ax.YTick = [0 .5 1 1.5 2 2.5];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% GLOBAL ERROR:LOG PLOT (BESSEL K)
clear all; clc;
x = linspace(-30,0,100); tt = 10.^x;
savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha1over5/Examples/';

load([savepath 'gauss_m0s1Xgauss_m0s1_alpha025']) %'W_vec', 'W_kj', 'k_start','k_end','j_start','j_end')
alpha = 0.25;
p025 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

load([savepath 'gauss_m0s1Xgauss_m0s1_alpha03']) %'W_vec', 'W_kj', 'k_start','k_end','j_start','j_end')
alpha = 0.3;
p03 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

load([savepath 'gauss_m0s1Xgauss_m0s1_alpha035']) %'W_vec', 'W_kj', 'k_start','k_end','j_start','j_end')
alpha = 0.35;
p035 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

load([savepath 'gauss_m0s1Xgauss_m0s1_alpha045']) %'W_vec', 'W_kj', 'k_start','k_end','j_start','j_end')
alpha = 0.45;
p045 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

%besselK
s_x = 1; s_y = 1;
p_trueK = @(x) besselk( 0, abs(x)/(s_x * s_y) ) / (pi * s_x * s_y);
pt = p_trueK(tt);


%color figure
figure; 
plot(x,log10(abs(p045-pt)./pt + 1e-20));hold on
plot(x,log10(abs(p035-pt)./pt + 1e-20));hold on
plot(x,log10(abs(p03-pt)./pt + 1e-20));hold on
plot(x,log10(abs(p025-pt)./pt + 1e-20));hold on

%black figure
figure; 
plot(x,log10(abs(p045-pt)./pt + 1e-20),'k-');hold on
plot(x,log10(abs(p035-pt)./pt + 1e-20),'k--');hold on
plot(x,log10(abs(p03-pt)./pt + 1e-20),'k-.');hold on
plot(x,log10(abs(p025-pt)./pt + 1e-20),'k:');hold on

axis([-30 0 -16 0])
% %setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-30:5:0]; ax.YTick = [-16:2:0];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% GLOBAL ERROR (BESSEL K)
tt = linspace(-10,10,5000);

p = eval_GMR(W_vec,tt,alpha);
% eval using support of basis functions
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

%p_near_zero = eval_GMR(W_vec,0.1)
%p2 = eval_GMR(W_vec2,t);

%besselK
s_x = 1; s_y = 1;
p_trueK = @(x) besselk( 0, abs(x)/(s_x * s_y) ) / (pi * s_x * s_y);

figure; plot(tt,p2-p_trueK(tt),'k');

%  axis([-10 10 -1.5e-15 2.5e-15])
% % %setting axis tick locations, labels and font size
% ax = gca; ax.XTick = [-10 -5 0 5 10]; ax.YTick = [-1 0 1 2]*1e-15;
% ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% LOCAL (near singularity ) ERROR (BESSEL K) 
tt = linspace(-2e-10, 2e-10,5000);


p = eval_GMR(W_vec,tt,alpha);

p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

%p_near_zero = eval_GMR(W_vec,0.1)
%p2 = eval_GMR(W_vec2,t);

%besselK
s_x = 1; s_y = 1;
p_trueK = @(x) besselk( 0, abs(x)/(s_x * s_y) ) / (pi * s_x * s_y);

%figure; plot(tt,(p-p_trueK(tt)),'k');
figure; plot(tt,(p2-p_trueK(tt)),'k');

%  axis([-2e-10 2e-10 -2.5e-14 2.5e-14])
% % %setting axis tick locations, labels and font size
% ax = gca; ax.XTick = [-2 -1 0 1 2]*1e-10; ax.YTick = [-2 -1 0 1 2]*1e-14;
% ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% SINGULARITY - 1 (where round tip is seen)
tt = linspace(-1e-28,1e-28,5000);

p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p2,'k');

axis([-1e-28 1e-28 20.5 21.9])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-1 -.5 0 .5 1]*1e-28; ax.YTick = [20.6 21.2 21.8];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% SINGULARITY - 2 (where sharp tip is seen)
tt = linspace(-1e-27,1e-27,5000);

p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p2,'k');

axis([-1e-27 1e-27 19.8 21.9])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-1 -.5 0 .5 1]*1e-27; ax.YTick = [20 20.6 21.2 21.8];
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


