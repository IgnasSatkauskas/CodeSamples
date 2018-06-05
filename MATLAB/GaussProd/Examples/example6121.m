%example 6226
clear all; clc;
%close all;

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

format long

alpha = 0.25;
%alpha = 1/3; %epsilon is 4.54 * 10^-10
%alpha = 1/5; %


% f(x) distribution can be arbitrary 
    
% f(x) is Gaussian
mu_x = 6;
s_x = 1;
f = @(x) ( 1 / (s_x*sqrt(2*pi)) ) * exp( -(x - mu_x).^2 / (2*s_x^2) );


% g(y) distribution is Gaussian 
mu_y = 2;
s_y = 1;


savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha13over30/Examples/';

%filename for saving runs

savename=['gauss_m6s1Xgauss_m2s1']

save_run = 0;

%COMPUTE w_kj 



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
%savename

%% Moments 
moment0 = 1
moment1 = mu_x * mu_y
moment2 = mu_x^2*s_y^2 + mu_y^2*s_x^2 + s_x^2*s_y^2 


moments(W_kj, j_start, j_end, k_start, k_end, alpha)


%PLOTTING-----------------------------------------------------------------

%% FACTORS
% f(x) is Gaussian
t = linspace(-10,20,10000);
figure; plot(t,f(t),'k'); 
axis([-10 20 0 .5])

% g(x) is Gaussian 
g = @(x) (1/(sqrt(2*pi)*s_y)) * exp(-(x-mu_y).^2 / (2*s_y^2));
hold on; plot(t, g(t),'k'); 

%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-10 -5 0 5 10 15 20]; ax.YTick = [0 .1 .2 .3 .4 .5];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% PRODUCT
tt = linspace(-20,40,5000);

p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p2,'k');
axis([-20 40 0 .1])

% %setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-20 -10 0 10 20 30 40]; ax.YTick = [0 .05 .1];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;




%% SINGULARITY 1
tt = linspace(-1e-6,1e-6,10000);

p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p2,'k');

 axis([-1e-6 1e-6 9.272957*1e-3 9.272969*1e-3])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-1 -.5 0 .5 1]*1e-6; ax.YTick = [9.272958 9.272963 9.272968]*1e-3;
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% SINGULARITY 2
tt = linspace(-1e-6,1e-6,10000);

p = eval_GMR(W_vec,tt,alpha);
p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p2,'k');

% axis([-2e-10 2e-10 9.272967*1e-3 9.272974*1e-3])
% %setting axis tick locations, labels and font size
% ax = gca; ax.XTick = [-2 -1 0 1 2]*1e-10; ax.YTick = [9.272968 9.272971 9.272974]*1e-3;
% ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

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

colormap gray

%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])

%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-100 -50 0 50 100]; ax.YTick = [0 30 60 90];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;


%% CHECK DIFFERENCE BETWEEN TWO REPRESENTATIONS  62 and 26
clear all; clc;
path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

%alpha=0.2;
alpha=13/30;

%load 
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/';
savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha13over30/Examples/';

savename1 = 'gauss_m6s1Xgauss_m2s1' ;
D1=load([savepath savename1]);
W_vec1=D1.W_vec;
W_kj1 = D1.W_kj;

%
k_start = D1.k_start; k_end = D1.k_end; j_start = D1.j_start; j_end = D1.j_end;


savename2 = 'gauss_m2s1Xgauss_m6s1' ;
D2=load([savepath savename2]);
W_vec2=D2.W_vec;
W_kj2 = D2.W_kj;


%evaluate product distributions
t = linspace(-20,50,10000);
t = linspace(-2e-6,2e-6,10000);

p62 = eval_GMR(W_kj1,t,alpha); 
p26 = eval_GMR(W_kj2,t,alpha);

figure; plot(t,p62-p26); title('difference bw 62 and 26')


