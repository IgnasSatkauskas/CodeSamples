%LaplaceXGumbel example
% Gumbel as f(x)
% Laplace as sum of exponentianls

clear all;
path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );
addpath( [path '/Examples'] );
format long

%Laplace as sum of exponentials
% in LyxNotes/LaplaceXGumbel.lyx

%first parameters of Laplace PDF:
mu = 0;
b = 1;

%figure; ; plot(tt,g_true(tt)); title('Laplace PDF')


% %dicretization step:
% h = 0.1;
% % discretization points t_n
% n = 150; %number of points on each side of mu
% % n points on each side of mu, h distance apart:
% t_n = [-n*h+mu:h:mu mu+h:h:mu+n*h];

t_start = -10;
t_end = 5;
n_points = 50;
h = (t_end - t_start)/(n_points-1)
t_n =  linspace(t_start,t_end,n_points);

%coeficients c_n
c_n = ( h/(4*b^2*sqrt(pi)) ) * ( exp(-exp(t_n)/(4*b^2) + t_n/2) );

%these go into computation of w_kj (have 1 / (signa*sqrt(2pi)) factored out)
cc_n = ( h/(4*b^2) ) * ( exp( -exp(t_n)/(4*b^2) + t_n ) );

% sigma_n's : sigmas of the Gaussians in Laplace expanssion
s_n = (1/sqrt(2)) * exp(t_n/2);

% resulting approximation of Laplace
g = @(x) sum_gauss(row(c_n), col(s_n), mu, x);

%------------------plotting--------------------------------------
tt = linspace(-1e-6,1e-6,1000);
t = linspace(-20,20,1000);
%Laplace true
g_true = @(x) (1/(2*b)) * (exp(-abs(x-mu)/b));
%plot Laplace true

% %plot approximation
% figure; plot(tt,g(tt)); title('Laplace Approximation')

% %check some point
% g_true(mu)
% g(mu)

%plot the error local
figure; plot(tt,g_true(tt)-g(tt)); title('Error local')
%plot the error global
figure; plot(t,g_true(t)-g(t)); title('Error global')

% total integral check
total_integral = sqrt(2*pi) * row(c_n) * col(s_n)

%first and last sigmas
sigma_1=s_n(1)
sigma_end=s_n(end)


%% ----------------------------------------------------------
% Now find PDF of Laplace X Gumbel 

%% 

alpha = 0.25;
% f(x) is Extreme Value Distribution
mu_x = 3; % center
b = 1; % shape
f = @(x) (1/b) * exp(-(x-mu_x)/b -exp(-(x-mu_x)/b));

%  W_kj size 
k_start = -100;
k_end = 100;
j_start = -5;
j_end = 100;

%W_kj = zeros(j_end-j_start+1,k_end-k_start+1);

% W_kj with snapshots
W_kj = zeros(j_end-j_start+1,k_end-k_start+1,length(c_n));

tic
% "multiply" Gumbel by each Gaussian in Laplace expansion
for n = 1:length(c_n)
    n
% g(x) is Gaussian 
mu_y = mu; %mean is the same for all Gaussians in expansion
s_y = s_n(n); % sigma's of Gaussians in expansion



tol_integral =10^-8;
%tol_weight = 10^-20;


c1 = 1 / (sqrt(2*alpha)*s_y) ; % constant repeats in integrand definition 4 times; precompute
j_i = 0; %W_kj row index
    
for j = j_start:j_end
    %j
    j_i = j_i+1;
    k_i = 0; % W_kj column index
    for k = k_start:k_end
        k_i = k_i+1;
        %integrand(tau,f,s_y,mu_y,j,k)
        integ = @(tau) integrand(tau,f,s_y,mu_y,alpha,j,k);
        %AdaptiveIntegrator
        %w_kj = ( ( 2^(-j/2)*log(2) ) / (sqrt(2*pi)*s_y) ) * AdaptiveIntegrator(integ,0,1,tol_integral);
        w_kj = ( ( 2^(-j/2)*log(2) ) * c_n(n) ) * AdaptiveIntegrator(integ,0,1,tol_integral);
     
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
        %W_kj(j_i,k_i) = W_kj(j_i,k_i) + cc_n(n) * w_kj;
        %W_kj(j_i,k_i) = W_kj(j_i,k_i) +  w_kj;
        W_kj(j_i,k_i,n) = W_kj(j_i,k_i,n) +  w_kj;
        
    end
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

%
% sum along 3rd dim
W_kj1 = sum(W_kj,3);


% total integral check
[j_i_end,k_i_end] = size(W_kj1);
total_integral =0;
for j_i = 1:j_i_end
    j = j_start + (j_i -1);
    total_integral = total_integral + sum(W_kj1(j_i,:))*2^(-j/2);
end
total_integral
%savename

%%
savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Examples/Runs/';
savename=['Laplace_m0b1XGumbel_mneg2b1_snapshots'];
save_run = 1;

if save_run
    %save([savepath savename], 'W_vec')
    save([savepath savename], 'W_vec', 'W_kj', 'k_start','k_end','j_start','j_end', 'alpha', 'W_kj1')
    savename
end





%% -----------------PLOTTING----------------------------

%% Load before ploting
loadpath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Examples/Runs/';
loadname=['Laplace_m0b1XGumbel_mneg1b1']
load([loadpath loadname])


%% PRODUCT
tt = linspace(-20,20,5000);

p2 = eval_GMR2(W_kj1,j_start,j_end,k_start,k_end,tt,alpha);
%p2 = eval_GMR2(W_kj(:,:,1),j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p2,'k');

% axis([-20 20 0 .5])
% % %setting axis tick locations, labels and font size
% ax = gca; ax.XTick = [-20 -10 0 10 20]; ax.YTick = [0 .1 .2 .3 .4 .5];
% ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;


%% COEFFICIENTS

figure;
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj1+1e-32),-32:.25:-2); %title('GMR matrix of coefficients (log scale)')
%[~,h]=contourf(X,Y,log10(W_kj(:,:,16)+1e-32),-32:.25:-2); %title('GMR matrix of coefficients (log scale)')
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






