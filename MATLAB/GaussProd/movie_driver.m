% SAVING MANY RUNS FOR MOVIES


%NEW FORMULATION (WITHOUT S-INTEGRAL)
clear all; clc;
%close all;

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

format long
alpha = 13/30;

% f(x) distribution can be arbitrary 
loops = 201;
for run_num = 1:loops

run_num
    
start_point = -5;
% % f(x) is Gaussian
% mu_x = 0 + (run_num-1)/5;
%  mu_x = start_point + (run_num-1)/10; % center
%  s_x = 0.3;
  mu_x = 5;
 s_x = 0.3;
 f = @(x) ( 1 / (s_x*sqrt(2*pi)) ) * exp( -(x - mu_x).^2 / (2*s_x^2) );

% % f(x) is a Laplace distribution
% %mu = 0 + (run_num-1)/5;
%   mu = start_point + (run_num-1)/10; % center
% b = 1;
% f = @(x) (1/(2*b)) * exp(-abs(x-mu)/b);

% % f(x) is Cauchy distribution
% %mu = 0 + (run_num-1)/5; %center
%  mu = start_point + (run_num-1)/10; % center
% b = 1; % shape parameter
% f = @(x) ( b^2 ./ ((x-mu).^2 + b^2) ) / (pi*b);

% % f(x) is Extreme Value Distribution
% mu = start_point + (run_num-1)/10; % center
% b = 1; % shape
% f = @(x) (1/b) * exp(-(x-mu)/b) .* exp( -exp(-(x-mu)/b) );

% % f(x) distribution as product of two other distributions 
% loadpath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha13over30/';
% loadname = 'gauss_m2s1Xlaplace_m1b1' %Gaussian and Laplace product
% load([loadpath loadname])
% W_vec_f = W_vec;
% clear W_vec;
% f = @(x) eval_GMR(W_vec_f,x,alpha);

% g(y) distribution is Gaussian 
mu_y = start_point + (run_num-1)/10; % center
 s_y = 0.3;
% mu_y = 5;
% s_y = 0.3;

% %flip order of multiplication
% mu_temp = mu_x; s_temp = s_x;
% mu_x = mu_y; s_x = s_y;
% mu_y = mu_temp; s_y = s_temp;



savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Movies/';

%savename=['gauss0Xgauss' num2str(run_num)]
%savename=['gauss0Xlaplace' num2str(run_num)]
%savename=['gauss2Xlaplace' num2str(run_num)]
%savename=['gauss0Xcauchy' num2str(run_num)] % cauchy needs j_star = -50
%savename=['gauss2Xcauchy' num2str(run_num)]
%savename=['gauss1Xcauchy' num2str(run_num)]
%savename=['gauss1s15Xcauchy' num2str(run_num)]
%savename=['gauss2Xextreme' num2str(run_num)] % extreme needt j_start = -5  
%savename=['gaussANDlaplaceXgauss' num2str(run_num)]
%savename=['gauss_m2s1Xextreme_b1/' 'frame' num2str(run_num)] % extreme needs j_start = -5 
%savename=['gauss_m2s1Xcauchy_b1/' 'frame' num2str(run_num)] % cauchy needs j_star = -50
%savename=['gauss_m2s1Xlaplace_b1/' 'frame' num2str(run_num)] %
%savename=['gauss_m2s1Xgauss_s1/' 'frame' num2str(run_num)] %
%savename=['gauss_m4s0p3Xgauss_s0p3/' 'frame' num2str(run_num)] %
savename=['gauss_s0p3Xgauss_m5s0p3/' 'frame' num2str(run_num)] %

save_run = 1;

%COMPUTE w_kj 
%alpha = 0.2;
%alpha = 1/3; %epsilon is 4.54 * 10^-10
alpha = 13/30; %epsilon ~10^-8

% compute w_kj
k_start = -200;
k_end = 200;
j_start = -10;
j_end = 100;

W_kj = zeros(j_end-j_start+1,k_end-k_start+1);



tol_integral =10^-14;
%tol_weight = 10^-20;


tic
k_i = 0; % W_kj column index
c1 = 1 / (sqrt(2*alpha)*s_y) ; % constant repeats in integrand definition 4 times; precompute

for k = k_start:k_end
    k_i = k_i+1;
    k;
    j_i = 0; %W_kj row index
    for j = j_start:j_end
     
        j_i = j_i+1;
        
        
        
        %integrand(tau,f,s_y,mu_y,j,k)
        integ = @(tau) integrand(tau,f,s_y,mu_y,alpha,j,k);

       
        %AdaptiveIntegrator
        w_kj = ( ( 2^(-j/2)*log(2) ) / (sqrt(2*pi)*s_y) ) * AdaptiveIntegrator(integ,0,1,tol_integral);
        
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
end


% total integral check
[j_i_end,k_i_end] = size(W_kj);
total_integral =0;
for j_i = 1:j_i_end
    j = j_start + (j_i -1);
    total_integral = total_integral + sum(W_kj(j_i,:))*2^(-j/2);
end
total_integral
savename

end % for run_num


%% Construct the movie frames and play it
alpha = 13/30;
savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Movies/';
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Movies/Png/gauss2Xextreme/';



F(loops) = struct('cdata',[],'colormap',[]);

for run_num = 1:loops
    %load run
    %savename=['gauss0Xgauss' num2str(run_num)];
    %savename=['gauss2Xlaplace' num2str(run_num)];
    %savename=['gauss0Xcauchy' num2str(run_num)];
    %savename=['gauss2Xcauchy' num2str(run_num)];
    %savename=['gauss1Xcauchy' num2str(run_num)];
    %savename=['gauss1s15Xcauchy' num2str(run_num)];
    %savename=['gaussANDlaplaceXgauss' num2str(run_num)];
    %savename=['gauss_m2s1Xextreme_b1/' 'frame' num2str(run_num)];
    %savename=['gauss_m2s1Xcauchy_b1/' 'frame' num2str(run_num)];
    %savename=['gauss_m2s1Xlaplace_b1/' 'frame' num2str(run_num)];
    %savename=['gauss_m4s0p3Xgauss_s0p3/' 'frame' num2str(run_num)];
    savename=['gauss_s0p3Xgauss_m5s0p3/' 'frame' num2str(run_num)];
    
    load([savepath savename])%  'W_vec', 'W_kj', 'k_start','k_end','j_start','j_end'
    tt = linspace(-40,100,20000);
    product_dist = eval_GMR(W_vec,tt,alpha);
    
    %plot factors
        % f - moving Gaussian
        %mu_x = 0 + (run_num-1)/5;
        mu_x = start_point + (run_num-1)/10; % center
        s_x = 0.3;
        f = @(x) ( 1 / (s_x*sqrt(2*pi)) ) * exp( -(x - mu_x).^2 / (2*s_x^2) );
    
%     % f(x) is a Laplace distribution
%     %mu = 0 + (run_num-1)/5;
%       mu = start_point + (run_num-1)/10; % center
%     b = 1;
%     f = @(x) (1/(2*b)) * exp(-abs(x-mu)/b);

% % f(x) is Cauchy distribution
% %mu = 0 + (run_num-1)/5; %center
%  mu = start_point + (run_num-1)/10; % center
% b = 1; % shape parameter
% f = @(x) ( b^2 ./ ((x-mu).^2 + b^2) ) / (pi*b);

% % f(x) is Extreme Value Distribution
% mu = start_point + (run_num-1)/10; % center
% b = 1; % shape
% f = @(x) (1/b) * exp(-(x-mu)/b) .* exp( -exp(-(x-mu)/b) );

% %f is product of other two distributions
% f = @(x) eval_GMR(W_vec_f,x,alpha);
    
    % g(y)  is fixed Gaussian
    mu_y = 5;
    s_y = 0.3;
    g = @(x) (1/(sqrt(2*pi)*s_y)) * exp(-(x-mu_y).^2 / (2*s_y^2));
    
    t2 = linspace(-15,25,10000);
    subplot(3,1,1)
    plot(t2,g(t2),t2,f(t2));
    
    subplot(3,1,2)
    plot(tt,product_dist);
    
    axis([-40 100 0 1])
    
    subplot(3,1,3)
    [X,Y]=meshgrid( k_start:k_end,j_start:j_end);
    % USE COUNTOURF INSTEAD OF SURF
    [~,h]=contourf(X,Y,log10(W_kj+eps),-16:.25:-2); %title('GMR matrix of coefficients (log scale)')
    set(h,'lineStyle','none');

    drawnow
    F(run_num) = getframe(gcf);
    
    filename = ['frame' num2str(run_num) '.png']
    saveas(gcf, [savepath 'gauss_s0p3Xgauss_m5s0p3/' filename])
    clear 'W_vec' 'W_kj' 'k_start' 'k_end' 'j_start' 'j_end'
    
end

%% play movie
fig = figure;
frames_per_sec = 10;
movie(fig, F,1, frames_per_sec)

%% save movie
save([savepath 'Fgauss2Xextreme'], 'F')

%% play saved movie
%load([savepath 'Fgauss1s15Xcauchy'])
%load([savepath 'Fgauss1Xcauchy'])
%load([savepath 'Fgauss2Xcauchy']) %similar to Fgauss1Xcauchy
%load([savepath 'Fgauss0Xcauchy'])
load([savepath 'Fgauss2Xlaplace'])
%load([savepath 'Fgauss2Xextreme'])


fig = figure;
frames_per_sec = 2;
movie(fig, F,1, frames_per_sec)

%% plot instance
    savename=['gauss1Xcauchy' num2str(6)];
    load([savepath savename])%  'W_vec', 'W_kj', 'k_start','k_end','j_start','j_end'
    tt = linspace(-10,10,10000);
    product_dist = eval_GMR(W_vec,tt,alpha);
    figure; plot(tt,product_dist); title(savename)
    %axis([-2 5 0 5])
    
 %% Plot coefficients
%old coefs
figure;
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj+eps),-16:.25:-2); %title('GMR matrix of coefficients (log scale)')
set(h,'lineStyle','none');

c=colorbar;
c.FontSize = 18;
c.Ticks=[-15:2:-3];

%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])   
title(savename)

