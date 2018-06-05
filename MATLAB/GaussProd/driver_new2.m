%NEW FORMULATION (WITHOUT S-INTEGRAL)
clear all; clc;
%close all;

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

format long

%alpha = 0.2;
%alpha = 1/3; %epsilon is 4.54 * 10^-10
alpha = 13/30; %epsilon


%f(x) distribution can be arbitrary 
    
% f(x) is Gaussian
mu_x = 2;
s_x = 1;
f = @(x) ( 1 / (s_x*sqrt(2*pi)) ) * exp( -(x - mu_x).^2 / (2*s_x^2) );

% % f(x) is a Laplace distribution
% mu = -2;
% b = 1;
% f = @(x) (1/(2*b)) * exp(-abs(x-mu)/b);





% % f(x) is Cauchy distribution
% mu = -1; %center
% b = 1; % scale parameter
% f = @(x) ( b^2 ./ ((x-mu).^2 + b^2) ) / (pi*b);

% % f(x) is Extreme Value Distribution
% mu = 0; % center
% b = 1; % shape
% f = @(x) (1/b) * exp(-(x-mu)/b) .* exp( -exp(-(x-mu)/b) );

% % f(x) distribution as product of two other distributions 
% loadpath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha13over30/';
% loadname = 'gauss_m2s1Xlaplace_m1b1' %Gaussian and Laplace product
% load([loadpath loadname])
% W_vec_f = W_vec;
% clear W_vec;
% f = @(x) eval_GMR(W_vec_f,x,alpha);



% % f(x) distribution as product of two 01 01 Gaussians
% load([path '/Runs/W0101_vec_14dig']);
% W0101_vec = W_vec;
% clear W_vec;
% f = @(x) eval_GMR(W0101_vec,x,alpha);

% g(y) distribution is Gaussian 
mu_y = 1;
s_y = 1;




%
%savepath to save runs
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/';
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha3/';
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha13over30/';
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Movies/';
savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha13over30/Examples/';

%filename for saving runs

%savename = ['M' strrep(num2str(mu_x),'.','p') 'S' strrep(num2str(s_x),'.','p') 'M' strrep(num2str(mu_y),'.','p') 'S' strrep(num2str(s_y),'.','p')];
%savename=['gauss_m2s1Xlaplace_m1b1']
%savename=['gauss_m2s1Xlaplace_m1b1Xgauss_m1s0p5']
%savename=['gauss_m2s1Xcauchy_mneg1b1']
savename=['gauss_m2s1Xlaplace_mneg2b1']

save_run = 0;

%COMPUTE w_kj 


% % define w_plus and w_minus
% w_plus = @(tau,k,j) f( (2^(2-j)*2.^(-tau)) / ( sqrt(2*alpha)*s_y) ) .* exp( ( (-alpha*4.^(tau-2)) ./ (1 - 4.^(tau-2)) ) .* (k - (mu_y * 2.^(-tau+2)) ./ (sqrt(2*alpha)*s_y) ).^2 );
% w_minus = @(tau,k,j) f( (-2^(2-j)*2.^(-tau)) / ( sqrt(2*alpha)*s_y) ) .* exp( ( (-alpha*4.^(tau-2)) ./ (1 - 4.^(tau-2)) ) .* (k + (mu_y * 2.^(-tau+2)) ./ (sqrt(2*alpha)*s_y) ).^2 );



% compute w_kj
k_start = -100;
k_end = 100;
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
    k
    j_i = 0; %W_kj row index
    for j = j_start:j_end
     
        j_i = j_i+1;
        
        
        %integrand = @(tau) ( w_plus(tau,k,j) + w_minus(tau,k,j) ) ./ ( sqrt(1-4.^(tau-2)) ) ;
        
        %w_plus = @(tau,k,j) f( (2^(2-j)*2.^(-tau)) / ( sqrt(2*alpha)*s_y) ) .* exp( ( (-alpha*4.^(tau-2)) ./ (1 - 4.^(tau-2)) ) .* (k - (mu_y * 2.^(-tau+2)) ./ (sqrt(2*alpha)*s_y) ).^2 );
        %w_minus = @(tau,k,j) f( (-2^(2-j)*2.^(-tau)) / ( sqrt(2*alpha)*s_y) ) .* exp( ( (-alpha*4.^(tau-2)) ./ (1 - 4.^(tau-2)) ) .* (k + (mu_y * 2.^(-tau+2)) ./ (sqrt(2*alpha)*s_y) ).^2 );
        
        %integrand = @(tau) ( f( c1*2^(2-j)*2.^(-tau) ) .* exp( ( (-alpha*4.^(tau-2)) ./ (1 - 4.^(tau-2)) ) .* (k - c1*mu_y * 2.^(-tau+2) ).^2 ) ...
        %                  + f( -c1*2^(2-j)*2.^(-tau) )  .* exp( ( (-alpha*4.^(tau-2)) ./ (1 - 4.^(tau-2)) ) .* (k + c1*mu_y * 2.^(-tau+2) ).^2 ) ) ./ ( sqrt(1-4.^(tau-2)) );
        
        %integrand(tau,f,s_y,mu_y,j,k)
        integ = @(tau) integrand(tau,f,s_y,mu_y,alpha,j,k);

        %quad
        %w_kj = ( ( 2^(-j/2)*log(2) ) / (sqrt(2*pi)*s_y) ) * quad(integrand,0,1,tol_integral);
        
        %quadgk
        %w_kj = ( ( 2^(-j/2)*log(2) ) / (sqrt(2*pi)*s_y) ) * quadgk(integrand,0,1);
        
        %AdaptiveIntegrator
        w_kj = ( ( 2^(-j/2)*log(2) ) / (sqrt(2*pi)*s_y) ) * AdaptiveIntegrator(integ,0,1,tol_integral);
        
        W_kj(j_i,k_i) = w_kj;
        
%         if abs(w_kj) > tol_weight
%             W_kj(j_i,k_i) = w_kj;
%         end
        
        
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



% save([path '/Runs/W010101_vec_14dig'],'W_vec')
% save([path '/Runs/W010101_14dig'],'W_kj', 'k_start', 'k_end', 'j_start', 'j_end')
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


moments(W_kj, j_start, j_end, k_start, k_end, alpha, 2)




%% Plotting


tt = linspace(-30,30,10000);

p = eval_GMR(W_vec,tt,alpha);
%p_near_zero = eval_GMR(W_vec,0.1)
%p2 = eval_GMR(W_vec2,t);

%PRODUCT 
figure; plot(tt,p,'k');

% %setting axis tick locations, labels and font size
% ax = gca; ax.XTick = [-1 -0.5 0 0.5 1]*1e-6; ax.YTick = [9.27296 9.27297]*1e-3;
% ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%title('Product Distribution')
%title(['Product distribution  ' savename]);
%figure; plot(t,p-p2); title(['error  ' savename]);

% plot f(x)

%% FACTORS
% f(x) is a Laplace distribution
mu = -2;
b = 1;
f = @(x) (1/(2*b)) * exp(-abs(x-mu)/b);
t = linspace(-20,20,10000);
figure; plot(t,f(t),'k'); 
axis([-20 20 0 .6])


% g(x) is Gaussian 
mu_y =2; s_y = 1;
g = @(x) (1/(sqrt(2*pi)*s_y)) * exp(-(x-mu_y).^2 / (2*s_y^2));
hold on; plot(t, g(t),'k'); 

%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-20 -10 0 10 20]; ax.YTick = [0 .2 .4 .6];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%%
%title('factors of the product')

%
% COEFFICIENTS
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
