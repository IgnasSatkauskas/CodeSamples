%% Movie from figs
clear all; clc;

s_y = 1;
mu_y = 1;
g = @(x,j) (1/(sqrt(2*pi)*s_y)) * exp(-(x-j*mu_y).^2 / (2*s_y^2));


loops = 20;
F(loops) = struct('cdata',[],'colormap',[]);
for j = 1:loops
    tt = linspace(-20,40,10000);
    plot(tt,g(tt,j));
    drawnow
    F(j) = getframe(gcf);
end

%% play movie
fig = figure;
movie(fig,F,2)


%% Load for plotting
clear all; %close all
format longe
path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
%path at cu machine
% path = '/home/satkausk/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

% load representation 
%loadpath = [path '/Runs/Alpha3/'];
%loadpath = [path '/Runs/Alpha3/Neg_coef/'];
%loadpath = [path '/Runs/Alpha13over30/'];
loadpath = [path '/Runs/Alpha1over4/'];
%loadname = 'M2S0p1M6S0p1no5no4';
%loadname = 'M2S0p1M6S0p1'
%loadname = 'M2S0p1M6S0p1no5'
%loadname = 'M0S1M0S1';
%loadname = 'going_up_from_2_to_4';
%loadname = 'going_up_from_2_to_5';
loadname = 'going_down_from_4_to_3';
%loadname = 'going_down_from_3_to_2';
%loadname = 'going_down_from_5_to_4';
%loadname = 'going_down_from_5_to_4_to_3';
load([loadpath loadname]); % 'W_vec', 'W_kj', 'k_start','k_end','j_start','j_end'



alpha = 1/4;


%% Plot coefficients
figure;
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(abs(W_kj)+eps),-16:.25:-2); %title('GMR matrix of coefficients (log scale)')
set(h,'lineStyle','none');

c=colorbar;
c.FontSize = 18;
c.Ticks=[-15:2:-3];



loadname4title = strrep(loadname, '_', ' ')
title(['Coefficients for  ' loadname4title ])
axis([k_start k_end j_start j_end])


%% Plot functions and errors

%basis function
phi = @(x,j,k) 2^(j/2) * sqrt(alpha/pi) * exp(-alpha*(2^j*x-k).^2);

t = linspace(-10,10,20000);



old = eval_GMR(W_vec,t,alpha);
true = phi(t,2,1);
error = old - true;

% error of new and old with actual basis funtion we startet with phi(x,2,1)
figure; plot(t,error); title(['Error bw phi and  ' loadname4title ]);








%% CHECK m0 and M00
clear all; 
%regular
alpha = 1/2;
Pi = pi;
k = -20:1:20; %m-knot's k
one = 1;
two = 2;
three = 3;


m0vpa = @(p)  sqrt(alpha/(3*Pi)) * exp(-(alpha*k.^2)/3) * ( exp(-2*Pi*1i*col(p)*k) ).' ;
M00vpa = @(p)  (  (sqrt(alpha/(3*Pi)) * exp(-(alpha*k.^2)/3) * ( exp(-2*Pi*1i*col(p)*k) ).' )  .* ( exp(-(alpha*k.^2)/2) * ( exp(2*Pi*1i*col(p)*k) ).' ) ) ./ (exp(-(alpha*k.^2)/2) * ( exp(4*Pi*1i*col(p)*k) ).' ) ;

% check if zero
pt = linspace(-0.75,0.75,10000);
zero_check = M00vpa(pt) .* m0vpa(pt) + M00vpa(pt+1/2) .* m0vpa(pt+1/2) - 1;

figure; plot(pt,real(zero_check)); title(num2str(alpha))
%figure; plot(pt,imag(zero_check)); title('imag')

%% VPA ----------------------------------------------------------------
clear all; 
%precision
d=32;
%alpha = vpa(0.2,d)
alpha = vpa(1/2,d);
Pi = vpa(pi,d);
k = -20:1:20; %m-knot's k
one = vpa(1,d);
two = vpa(2,d);
three = vpa(3,d);
four = vpa(4,d);

m0vpa = @(p) double( vpa( sqrt(alpha/(three*Pi)) * exp(-(alpha*k.^two)/three) * ( exp(-two*Pi*1i*col(p)*k) ).' ,d) );
M00vpa = @(p) double( vpa( (  (sqrt(alpha/(three*Pi)) * exp(-(alpha*k.^two)/three) * ( exp(-two*Pi*1i*col(p)*k) ).' )  .* ( exp(-(alpha*k.^two)/two) * ( exp(two*Pi*1i*col(p)*k) ).' ) ) ./ (exp(-(alpha*k.^two)/two) * ( exp(four*Pi*1i*col(p)*k) ).' ) ,d) );

% check if zero
pt = linspace(-0.75,0.75,1000);
zero_check = M00vpa(pt) .* m0vpa(pt) + M00vpa(pt+1/2) .* m0vpa(pt+1/2) - 1;

figure; plot(pt,real(zero_check));title('VPA: 0.5 ')
%figure; plot(pt,imag(zero_check),'ko');title('imag')

%% GOING UP

beta = 1;
s = 3;

f1 = @(x) exp(-beta*(x-s).^2);

%plot
t = linspace(-5,9,1000);
figure; plot(t,f(t))

phi = @(x,j,k) 2^(j/2) * sqrt(alpha/pi) * exp(-alpha*(2^j*x-k));

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
samples = linspace(-t_point,t_point,N_samples).';

%right hand side - evaluate function at sample points
rhs = eval_GMR(W_vec_scales, samples, alpha) / (sqrt(alpha/pi));


N_coefs_new = 70; % number of new coefs is 2N+1
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
delta = 10^-6; %away from zero
t = [ linspace(-80,-delta,20000) linspace(delta, 80,20000)];
f = eval_GMR(W_vec_scales,t,alpha);
figure; plot(t,f); title('f')

f_new = eval_GMR(W_vec_0scale_new,t,alpha);
figure; plot(t,f_new); title('f new')

figure; plot(t,f - f_new); title('error')

f_end = f(end)







