%% REDUCTION: GETTING BACK THE COEFS WITH LINEAR SOLVE 
% RANDOM COEFFICIENTS TEST CASE 

format longe
clear all; clc

%number of coefs is 2N+1
N = 30;

%pick random coefficients
coefs = rand(2*N+1,1);

%% make gaussian envelope for coefs
sigma =5;
xx = (-N:N) ./ sigma;
factor = exp(-xx.^2);
coefs = factor.' .* coefs;

%%
alpha = 13/30;
shifts = -N:N;

%function  temp(x, coefs, shifts)
%sum of Gaussians with shifts and coefs

x = linspace(-N-10,N+10,10000);
figure(1); plot(coefs,'ko')
figure(2); plot(x,temp(x,coefs,shifts)); title('original: sum of few gaussians')

%check
last_value = temp(N+10,coefs,shifts)
last_coef = coefs(end)

%%
%sample points
N_samples = 5000;
samples = linspace(-N-10,N+10,N_samples);

%right hand side
rhs = temp(samples,coefs,shifts).';

% construct a matrix
A = zeros(length(samples), length(coefs));

for k = 1:length(coefs)
    
    A(:,k) = exp(-alpha*(samples-shifts(k)).^2);
    
end

coefs_new = A\rhs

difference = coefs_new-coefs

[coefs coefs_new]

max_error = max(abs(coefs_new-coefs))

%% lsqlin solve overdetermined system with constraints
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub)
options = optimoptions(@lsqlin,'Algorithm','active-set')
%x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options) 
%minimizes with an initial point x0 and the optimization options specified in options. 
%Use optimoptions to set these options. If you do not want to include an initial point, set the x0 argument to [].
x0 = [];
%%
Aeq = A; beq = rhs;
% lower bound of solution
lb = zeros(size(coefs));
%upperbound of solution
ub = Inf*ones(size(coefs));
coefs2 = lsqlin(A, rhs, [], [], Aeq, beq, lb, ub, x0, options);
[coefs coefs2]

%% REDUCTION: GETTING BACK THE COEFS WITH LINEAR SOLVE 
% TRUE COEFFICIENTS ON ZERO SCALE

% load representations
format longe
clear all; clc;
path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

%alpha=0.2;
alpha=13/30;

%load 
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/';
savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha13over30/';

%savename1 = 'M6S0p1M2S0p1' ;
savename1 = 'gauss_m2s1Xlaplace_m1b1Xgauss_m0s1.mat';
D1=load([savepath savename1]);
W_vec1=D1.W_vec;
W_kj1 = D1.W_kj;

%
k_start = D1.k_start; k_end = D1.k_end; j_start = D1.j_start; j_end = D1.j_end;

%savename2 = 'M2S0p1M6S0p1' ;
savename2 = 'Product3/gauss_m2s1Xlaplace_m1b1Xgauss_m0s1.mat';
D2=load([savepath savename2]);
W_vec2=D2.W_vec;
W_kj2 = D2.W_kj;

% % take one scale only - W_kj2 on zero scale
% W_kj = W_kj2; % which representation
% scale = 0;
% W_kj_1scale = zeros(size(W_kj1));
% W_kj_1scale(scale-j_start+1,:) =  W_kj(scale-j_start+1,:);
% % W_kj to W_vec
% W_vec_1scale = Wkj2Wvec(W_kj_1scale,k_start, j_start);


% take few scales only
W_kj = W_kj2; % which representation
scale_start = 0;
scale_end = 4;
W_kj_scales = zeros(size(W_kj));
W_kj_scales(scale_start-j_start+1:scale_end-j_start+1,:) =  W_kj(scale_start-j_start+1:scale_end-j_start+1,:);
% W_kj to W_vec
W_vec_scales = Wkj2Wvec(W_kj_scales,k_start, j_start);



%% PLOTTING
% plot coefficients
figure
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);
[~,h]=contourf(X,Y,log10(W_kj_scales+eps),-16:.25:-2); title('Wkj2')
set(h,'lineStyle','none');
c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];
axis([k_start k_end j_start j_end])

%plot the function certain scales only
delta = 10^-6; %away from zero
t = [ linspace(-100,-delta,10000) linspace(delta, 100,10000)];
p = eval_GMR(W_vec_scales,t,alpha);
figure; plot(t,p)






%% LINEAR REDUCTION of one scale
format longe
%number of coefs is 2N+1
N = 133;
shifts = -N:N;
% take only coefs with shifts -N:N
W_vec_1scale = W_vec_1scale( abs(W_vec_1scale(:,2)) <= N , :);

coefs = W_vec_1scale(:,3);


%
%sample points
N_samples = 30000;
samples = linspace(-N-20,N+20,N_samples).';

%right hand side - evaluate function at sample points
rhs = eval_GMR(W_vec_1scale, samples, alpha) / (sqrt(alpha/pi));


N_coefs_new = 70; % want this  smaller that N
shifts_new = -N_coefs_new : N_coefs_new;

% construct a matrix
A = zeros(length(samples), length(shifts_new));

for k = 1:length(shifts_new)
    
    %make columns 
    A(:,k) =  exp(-alpha*(samples-shifts_new(k)).^2);
    
end

%% backslash
coefs_new = A\rhs

%% lsqlin solve overdetermined system with constraints
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub)
options = optimoptions(@lsqlin)%, 'TolFun', 10^-20)
%x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options) 
%minimizes with an initial point x0 and the optimization options specified in options. 
%Use optimoptions to set these options. If you do not want to include an initial point, set the x0 argument to [].
x0 = [];
%
%Aeq = A; beq = rhs;
Aeq = []; beq = [];
% lower bound of solution
lb = zeros(size(shifts_new));
%upperbound of solution
ub = Inf*ones(size(shifts_new));
coefs2 = lsqlin(A, rhs, [], [], Aeq, beq, lb, ub);
%

coefs_new = coefs2;

%% PLOTTING
    

% % plot coefs
% % 
% figure; plot(coefs_new); title('coefs new')
% figure; plot(coefs); title('coefs old')

% new W_vec representation
scale = 0;
W_vec_1scale_new = zeros(length(shifts_new),3);
W_vec_1scale_new(:,1) = scale * ones(length(shifts_new),1);
W_vec_1scale_new(:,2) = shifts_new;
W_vec_1scale_new(:,3) = coefs_new;

%
% plot old, new, and error
delta = 10^-6; %away from zero
t = [ linspace(-80,-delta,20000) linspace(delta, 80,20000)];
f = eval_GMR(W_vec_1scale,t,alpha);
figure; plot(t,f); title('f')

f_new = eval_GMR(W_vec_1scale_new,t,alpha);
figure; plot(t,f_new); title('f new')

figure; plot(t,f - f_new); title('error after reducing coefficients from 267 to 141')

f_end = f(end)


%% 
% basis elements comprising it 

phi = @(x,j,k) 2^(j/2) * sqrt(alpha/pi) * exp(-alpha*(2^j*x-k).^2);

W_vec_plot = W_vec_1scale_new(W_vec_1scale_new(:,3)>10^-8,:); %take only large coeficients for visualizing basis elements
tt=linspace(-40,40,10000);
figure
for plot_n = 1:length(W_vec_plot(:,1))
    basis_function = W_vec_plot(plot_n,3)* phi(tt, W_vec_plot(plot_n,1), W_vec_plot(plot_n,2));
    plot(tt,basis_function)
    hold on
end
% plot the original on top
hold on 
plot(tt, eval_GMR(W_vec_1scale,tt,alpha), 'k', 'LineWidth', 1)



%% LINEAR REDUCTION using few scales for rhs

% load representations
format longe
clear all; clc;
path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

alpha=13/30;

%load 
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/';
savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha13over30/';

%savename1 = 'M6S0p1M2S0p1' ;
savename1 = 'gauss_m2s1Xlaplace_m1b1Xgauss_m0s1.mat';
D1=load([savepath savename1]);
W_vec1=D1.W_vec;
W_kj1 = D1.W_kj;

k_start = D1.k_start; k_end = D1.k_end; j_start = D1.j_start; j_end = D1.j_end;

%savename2 = 'M2S0p1M6S0p1' ;
savename2 = 'Product3/gauss_m2s1Xlaplace_m1b1Xgauss_m0s1.mat';
D2=load([savepath savename2]);
W_vec2=D2.W_vec;
W_kj2 = D2.W_kj;

% drop scales below 0 scale for W_kj2
W_kj = W_kj2; % which representation
scale_start = 0; 
scale_end = 100;
W_kj2 = zeros(size(W_kj));
W_kj2(scale_start-j_start+1:scale_end-j_start+1,:) =  W_kj(scale_start-j_start+1:scale_end-j_start+1,:);
% W_kj to W_vec
W_vec2 = Wkj2Wvec(W_kj,k_start, j_start);



% take few scales only
W_kj = W_kj2; % which representation
scale_start = 0; 
scale_end = 10;
W_kj_scales = zeros(size(W_kj));
W_kj_scales(scale_start-j_start+1:scale_end-j_start+1,:) =  W_kj(scale_start-j_start+1:scale_end-j_start+1,:);
% W_kj to W_vec
W_vec_scales = Wkj2Wvec(W_kj_scales,k_start, j_start);



%% PLOTTING
% plot coefficients all scales 
figure
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);
[~,h]=contourf(X,Y,log10(W_kj2+eps),-16:.25:-2); title('Wkj2')
set(h,'lineStyle','none');
c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];
axis([k_start k_end j_start j_end])

% plot coefficients on some scales 
figure
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);
[~,h]=contourf(X,Y,log10(W_kj_scales+eps),-16:.25:-2); title('Wkj2')
set(h,'lineStyle','none');
c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];
axis([k_start k_end j_start j_end])

%plot the function - certain scales only
delta = 10^-6; %away from zero
t = [ linspace(-100,-delta,10000) linspace(delta, 100,10000)];
p = eval_GMR(W_vec_scales,t,alpha);
figure; plot(t,p)


%% RIGHT HAND SIDE and MATRIX for projecting on zero scale

%find t where function is negligible (for sampling)
for t_point = 0:400
    value = eval_GMR(W_vec_scales,t_point,alpha);
    if ( value < 10^-16 )
        break
    end
    t_point
end

%sample points
N_samples = 30000;
%samples = linspace(-t_point,t_point,N_samples).';
samples = linspace(-100,100,N_samples).';

%right hand side - evaluate function at sample points
rhs = eval_GMR(W_vec_scales, samples, alpha) / (sqrt(alpha/pi));


N_coefs_new = 80; % number of new coefs is 2N+1
shifts_new = -N_coefs_new : N_coefs_new;

% construct the matrix
A = zeros(length(samples), length(shifts_new));

for k = 1:length(shifts_new)
    
    %make columns 
    A(:,k) =  exp(-alpha*(samples-shifts_new(k)).^2);
    
end

%% backslash
coefs_new = A\rhs

%% lsqlin solve overdetermined system with constraints
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub)
options = optimoptions(@lsqlin)%, 'TolFun', 10^-20)
%x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options) 
%minimizes with an initial point x0 and the optimization options specified in options. 
%Use optimoptions to set these options. If you do not want to include an initial point, set the x0 argument to [].
x0 = [];
%
%Aeq = A; beq = rhs;
Aeq = []; beq = [];
% lower bound of solution
lb = zeros(size(shifts_new));
%upperbound of solution
ub = Inf*ones(size(shifts_new));
coefs2 = lsqlin(A, rhs, [], [], Aeq, beq, lb, ub);
%

coefs_new = coefs2;

%% PLOTTING
    
% new W_vec representation
scale = 0;
W_vec_0scale = zeros(length(shifts_new),3);
W_vec_0scale_new(:,1) = scale * ones(length(shifts_new),1);
W_vec_0scale_new(:,2) = shifts_new;
W_vec_0scale_new(:,3) = coefs_new;

%
% plot old, new, and error
delta = .1; %away from zero
t = [ linspace(-100,-delta,20000) linspace(delta, 100,20000)];
f = eval_GMR(W_vec_scales,t,alpha);
figure; plot(t,f); title('f')

f_new = eval_GMR(W_vec_0scale_new,t,alpha);
figure; plot(t,f_new); title('f new')

figure; plot(t,f - f_new); title('error')

f_end = f(end)


%% CHECK DIFFERENCE BW TWO REPRESENTATIONS : error when dropping below 0 scales

% load representations
format longe
clear all; clc;
path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

alpha=13/30;

%load 
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/';
savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha13over30/';

%savename1 = 'M6S0p1M2S0p1' ;
savename1 = 'gauss_m2s1Xlaplace_m1b1Xgauss_m0s1.mat';
D1=load([savepath savename1]);
W_vec1=D1.W_vec;
W_kj1 = D1.W_kj;

k_start = D1.k_start; k_end = D1.k_end; j_start = D1.j_start; j_end = D1.j_end;

%savename2 = 'M2S0p1M6S0p1' ;
savename2 = 'Product3/gauss_m2s1Xlaplace_m1b1Xgauss_m0s1.mat';
D2=load([savepath savename2]);
W_vec2=D2.W_vec;
W_kj2 = D2.W_kj;


% take 0:end scales
W_kj = W_kj2; % which representation
scale_start = 0; 
scale_end = 100;
W_kj_scales = zeros(size(W_kj));
W_kj_scales(scale_start-j_start+1:scale_end-j_start+1,:) =  W_kj(scale_start-j_start+1:scale_end-j_start+1,:);
% W_kj to W_vec
W_vec_scales = Wkj2Wvec(W_kj_scales,k_start, j_start);


% plot coefficients on some scales 
figure
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);
[~,h]=contourf(X,Y,log10(W_kj_scales+eps),-16:.25:-2); title('Wkj2')
set(h,'lineStyle','none');
c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];
axis([k_start k_end j_start j_end])

% plot coefficients on all scales 
figure
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);
[~,h]=contourf(X,Y,log10(W_kj2+eps),-16:.25:-2); title('Wkj2')
set(h,'lineStyle','none');
c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];
axis([k_start k_end j_start j_end])


%plot the error 
t = linspace(-200,200,20000);
p = eval_GMR(W_vec_scales,t,alpha); 
pp = eval_GMR(W_vec2,t,alpha); 
figure; plot(t,p-pp)

% %plot the function with dropped lower scales
% figure; plot(t,p)

