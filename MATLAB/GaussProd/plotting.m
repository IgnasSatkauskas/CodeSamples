%% PLOT USING savename
clear all; clc;
path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );
format longe

alpha =13/30;
%load 
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/';
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha3/';
%savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/';
savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/Alpha13over30/Examples/';

%savename = 'M8S1M12S1' ;
%savename = 'M8S1M14S1' ;
%savename = 'M6S1M2S1' ;
%savename = 'M2S0p1M6S0p1' ;
%savename = 'M6S0p1M2S0p1' ;
%savename = 'M14S1M8S1' ;
%savename = 'M2S1M1S1' ;
%savename = 'M0S1M0S1' ;
%savename = 'Alpha13over30/gauss_m2s1Xlaplace_m1b1';
%savename = 'Alpha13over30/gauss_m2s1Xlaplace_m1b1Xlaplace_m1b1.mat';
%savename = 'Alpha13over30/Product3/gauss_m2s1Xcauchy_m2b1Xlaplace_m1b1.mat'
%savename = 'Alpha13over30/Product3/gauss_m2s1Xlaplace_m1b1Xlaplace_m1b1.mat'
%savename = 'Alpha13over30/Product3/gauss_m2s1Xlaplace_m1b1Xgauss_m0s1.mat'
%savename = 'Alpha13over30/Product3/gauss_m2s1Xlaplace_m1b1Xgauss_m4s0p5.mat'
%savename = 'Alpha13over30/Reduction/test2_reduced_Prod3_gauss_m2s1Xlaplace_m1b1Xgauss_m0s1.mat'
%savename = 'Alpha13over30/gauss_m2s1Xlaplace_m1b1Xgauss_m0s1.mat'
%savename = 'Movies/gauss_m5s0p3Xgauss_s0p3/frame51.mat'
%savename = 'gauss_m2s1Xcauchy_mneg1b1';
savename = 'gauss_m2s1Xlaplace_mneg2b1';

% D2=load([savepath savename2]);
% W_vec2=D2.W_vec;

load([savepath savename]);
%'W_kj', 'k_start', 'k_end', 'j_start', 'j_end'
W_vec = W_vec(W_vec(:,3)>10^-16,:); %take only large coeficients 

%take real parts of W_kj ??
%W_kj = real(W_kj);

%plot the distribution
t = linspace(-30,30,20000);
%t = linspace(-1e-1,1e-1,10000); %for singularity at zero
%axis([-20 45 0 0.07])

p = eval_GMR(W_vec,t,alpha);
%p_near_zero = eval_GMR(W_vec,0.1)
%p2 = eval_GMR(W_vec2,t);

figure; plot(t,p,'k'); %title(['Product distribution  ' savename]);
%figure; plot(t,p-p2); title(['error  ' savename]);
axis([-20 20 0 .3])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-20 -10 0 10 20]; ax.YTick = [0 0.1 0.2 0.3];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

% total integral check
[j_i_end,k_i_end] = size(W_kj);
total_integral =0;
for j_i = 1:j_i_end
    j = j_start + (j_i -1);
    total_integral = total_integral + sum(W_kj(j_i,:))*2^(-j/2);
end
total_integral


%% PLOT COEFFICIENTS
figure

%surf(GMR); title('GMR matrix of coefficients')
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% surf(X,Y,log10(W_kj+eps)); %title('GMR matrix of coefficients (log scale)')
% %xlabel('translations')
% %ylabel('scales')
% %zlim([10^-10 10^-4])
% view([0 0 1])
% shading interp
% lighting phong

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(abs(W_kj)+eps),-16:.25:-2); %title('GMR matrix of coefficients (log scale)')
set(h,'lineStyle','none');

c=colorbar;
c.FontSize = 18;
c.Ticks=[-15:2:-3];

%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])

%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-100 -50 0 50 100]; ax.YTick = [-10 0 30 60 90];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

% % plot tips of basis functions
% centers = 2.^-(W_vec(:,1)) .* W_vec(:,2);
% heights = sqrt(2).^(W_vec(:,1)) * sqrt(1/(5*pi)) .* W_vec(:,3);
% figure
% plot(centers,heights,'ko'); title('basis functions (just the tips)')
% %layover the distribution
% hold on; 
% plot(t,p)


%% CHECK DIFFERENCE BETWEEN TWO REPRESENTATIONS  for linear reduction 
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

%% take one scale only - W_kj2 on zero scale
W_kj = W_kj2; % which representation
scale = 0;
W_kj_1scale = zeros(size(W_kj1));
W_kj_1scale(scale-j_start+1,:) =  W_kj(scale-j_start+1,:);
% W_kj to W_vec
W_vec_1scale = Wkj2Wvec(W_kj_1scale,k_start, j_start);

% plot coefficients
figure
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);
[~,h]=contourf(X,Y,log10(W_kj_1scale+eps),-16:.25:-2); title('Wkj2')
set(h,'lineStyle','none');
c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];
axis([k_start k_end j_start j_end])

%plot the function (one scale)
delta = 10^-6; %away from zero
t = [ linspace(-80,-delta,10000) linspace(delta, 80,10000)];
p = eval_GMR(W_vec_1scale,t,alpha); 
figure; plot(t,p, 'k'); title('Wvec2 only scale zero')


p_end = p(end)

% basis elements comprising it 

phi = @(x,j,k) 2^(j/2) * sqrt(alpha/pi) * exp(-alpha*(2^j*x-k).^2);

W_vec_plot = W_vec_1scale(W_vec_1scale(:,3)>10^-8,:); %take only large coeficients for visualizing basis elements
tt=linspace(-80,80,10000);
figure
for plot_n = 1:length(W_vec_plot(:,1))
    basis_function = W_vec_plot(plot_n,3)* phi(tt, W_vec_plot(plot_n,1), W_vec_plot(plot_n,2));
    plot(tt,basis_function)
    hold on
end
% plot the original on top
hold on 
plot(tt, eval_GMR(W_vec_1scale,tt,alpha), 'k', 'LineWidth', 1)
title('zero scale and basis elements comprising it')


%% plot representation on few scales W_vec1
figure
W_kj = W_kj1;
for scale = 0:3
    
    %take one scale
    W_kj_1scale = zeros(size(W_kj));
    W_kj_1scale(scale-j_start+1,:) =  W_kj(scale-j_start+1,:);
    % W_kj to W_vec
    W_vec_1scale = Wkj2Wvec(W_kj_1scale,k_start, j_start);
    
    delta = 10^-6; %away from zero
    t = [ linspace(-40,-delta,10000) linspace(delta, 40,10000)];
    p = eval_GMR(W_vec_1scale,t,alpha);
    
    
     plot(t,p); hold on
    
    
end

% draw main one on top
plot(t,eval_GMR(W_vec1, t,alpha), 'k', 'LineWidth',2); axis([-40,40,0,.1])
title('Wvec1 and what it is made of')

%% plot representation on few scales W_vec2
figure
W_kj = W_kj2;
for scale = 0:3
    
    %take one scale
    W_kj_1scale = zeros(size(W_kj));
    W_kj_1scale(scale-j_start+1,:) =  W_kj(scale-j_start+1,:);
    % W_kj to W_vec
    W_vec_1scale = Wkj2Wvec(W_kj_1scale,k_start, j_start);
    
    delta = 10^-6; %away from zero

    t = [ linspace(-40,-delta,10000) linspace(delta, 40,10000)];
    p = eval_GMR(W_vec_1scale,t,alpha);
    
    
     plot(t,p); hold on
    
    
end

% draw main one on top
plot(t,eval_GMR(W_vec1, t,alpha),'k', 'LineWidth',2); axis([-40,40,0,.1])
title('Wvec2 and what it is made of')




%%   zero out coefs above some scale
max_scale = 0;
W_kj1_one_scale = W_kj1;
W_kj2_one_scale = W_kj2;
W_kj1_one_scale(max_scale-j_start+1:end,:) = 0;
W_kj2_one_scale(max_scale-j_start+1:end,:) = 0;

% % just take max_scale i.e. zero out coefs below max_scale
W_kj1_one_scale(1:max_scale-j_start-1,:) = 0;
W_kj2_one_scale(1:max_scale-j_start-1,:) = 0;


% Make vector representation W_kj to W_vec (if dropping some coefs)
%lin_ind = find(W_kj_new > tol_coef);
lin_ind = find(abs(W_kj1_one_scale) > eps);
[j_i, k_i] = ind2sub(size(W_kj1_one_scale),lin_ind);

W_vec1_one_scale = zeros(length(lin_ind),3);
W_vec1_one_scale(:,1) = j_start + (j_i-1);
W_vec1_one_scale(:,2) = k_start + (k_i -1);
W_vec1_one_scale(:,3) = W_kj1_one_scale(lin_ind);

%lin_ind = find(W_kj_new > tol_coef);
lin_ind = find(abs(W_kj2_one_scale) > eps);
[j_i, k_i] = ind2sub(size(W_kj2_one_scale),lin_ind);

W_vec2_one_scale = zeros(length(lin_ind),3);
W_vec2_one_scale(:,1) = j_start + (j_i-1);
W_vec2_one_scale(:,2) = k_start + (k_i -1);
W_vec2_one_scale(:,3) = W_kj2_one_scale(lin_ind);






%% plotting the difference bw two representation : all coefs


delta = 10^-6; %away from zero

t = [ linspace(-20,-delta,10000) linspace(delta, 20,10000)];

%t = linspace(-1e-6,1e-6,10000); %for singularity at zero
%axis([-20 45 0 0.07])

p1 = eval_GMR(W_vec1,t,alpha);
p2 = eval_GMR(W_vec2,t,alpha);
figure; plot(t,p1); title('W-vec1')
figure; plot(t,p2); title('W-vec2')
figure; plot(t,p1-p2); title('difference bw two representations')
%title(['diff between ' savename1 ' and ' savename2])

%% PLOT COEFFICIENTS all coefs
figure
%surf(GMR); title('GMR matrix of coefficients')
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj1+eps),-16:.25:-2); title('Wkj1')
set(h,'lineStyle','none');

c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];

%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])

figure
%surf(GMR); title('GMR matrix of coefficients')
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj2+eps),-16:.25:-2); title('Wkj2')
set(h,'lineStyle','none');

c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];

%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])

%% PLOT COEFFICIENTS one scale
figure
%surf(GMR); title('GMR matrix of coefficients')
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj1_one_scale+eps),-16:.25:-2); title('Wkj1')
set(h,'lineStyle','none');

c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];

%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])

figure
%surf(GMR); title('GMR matrix of coefficients')
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj2_one_scale+eps),-16:.25:-2); title('Wkj2')
set(h,'lineStyle','none');

c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];

%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])

% smarter way
figure
%surf(GMR); title('GMR matrix of coefficients')
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj1_1scale+eps),-16:.25:-2); title('Wkj1-1scale')
set(h,'lineStyle','none');

c=colorbar; c.FontSize = 18; c.Ticks=[-15:2:-3];

%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])


%% PLOT PARTICULAR SCALES 
% make sure to run "zero out above and below" before
delta = 10^-6; %away from zero

t = [ linspace(-100,-delta,10000) linspace(delta, 100,10000)];


p1 = eval_GMR(W_vec1,t,alpha);
p2 = eval_GMR(W_vec2,t,alpha);
p3 = eval_GMR(W_vec1_1scale,t,alpha);

figure; plot(t,p1); title(['W-vec1 scale ' num2str(max_scale)])
figure; plot(t,p2); title(['W-vec2 scale ' num2str(max_scale)] )
figure; plot(t,p3); title(['W-vec1 scale ' num2str(scale)] )
figure; plot(t,p1-p2); title('difference bw two representations')

% plot  




%% RESOLVE SINGULARITY AT ZERO
clear all; clc;

%load 
savepath = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New/Runs/';
savename = 'M2S1M1S1' ;

load([savepath savename]);
%'W_kj', 'k_start', 'k_end', 'j_start', 'j_end'

%plot the distribution
%t = linspace(-20,20,20000);
a = -1.0e-25;
b = -a;

t = linspace(a,b,20000); %for singularity at zero

p = eval_GMR(W_vec,t);
%p2 = eval_GMR(W_vec2,t);

figure; plot(t,p,'k'); %title(['Singularity at 0 ' savename]);
%figure; plot(t,p-p2); title(['error  ' savename]);
%axis([a b 1.42 1.6])

%setting axis tick locations, labels and font size
% ax = gca; ax.XTick = [-0.8 -0.4 0 0.4 0.8]*1e-20; ax.YTick = [1.42 1.46 1.5 1.54 1.58];
ax = gca; ax.XTick = [-0.8 -0.4 0 0.4 0.8]*1e-25; ax.YTick = [1.5560268 1.5560272 1.5560276];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;



%%
%LOAD DISTRIBUTIONS

load([path '/Runs/W0101_vec_14dig']); % product of 2 Gaussians: W_vec
W2_vec = W_vec;
clear W_vec;

load([path '/Runs/W010101_vec_14dig']); % product of 3 Gaussians: W_vec
W3_vec = W_vec;
clear W_vec;

p = @(x) eval_GMR(W_vec,x);

%% PRODUCT OF 4 GAUSSIANS: compare with meijerG function
clear all; clc

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

load([path '/Runs/W01010101_vec_14dig']); % product of 4 Gaussians: W_vec

%p3 = @(x) eval_GMR(W_vec,x);

t = linspace(-4,4,1000);
p4 = eval_GMR(W_vec,t);

figure; plot(t,p4); title('our product of 4 Gaussians')

%meijerG (assuming sigma_x = sigma_y = sigma_z = sigma_w = 1)
n = 4;
z = t.^2 ./ 2^n;
const = 2^(n/2)*pi^(n/2);
p4_true = meijerGnum(n,z) / const;

figure; plot(t,p4_true); title('meijerG product of 4')

figure; plot(t,p4_true - p4); title('Error: product of 4')

%% PLOT COEFFICIENTS of 4 GAUSSIANS

clear all; clc

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

load([path '/Runs/W01010101_14dig']); % product of 4 Gaussians: W_kj

    figure
    %surf(GMR); title('GMR matrix of coefficients')
    [X,Y]=meshgrid( k_start:k_end,j_start:j_end);
    surf(X,Y,log10(W_kj+eps)); %title('GMR matrix of coefficients (log scale)')
    %xlabel('translations')
    %ylabel('scales')
    %zlim([10^-10 10^-4])
    view([0 0 1])
    shading interp
    lighting phong
    colorbar
    title('product of 4 Gaussians')
    axis([-100 100 -10 130])


%% PRODUCT OF 3 GAUSSIANS: compare with meijerG function
clear all; clc; close all

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

load([path '/Runs/W010101_vec_14dig']); % product of 3 Gaussians: W_vec

%p3 = @(x) eval_GMR(W_vec,x);

t = linspace(-10,10,1000);
p3 = eval_GMR(W_vec,t);

figure; plot(t,p3,'k'); %title('our product of 3 Gaussians')

%meijerG (assuming sigma_x = sigma_y = sigma_z = 1)
n = 3;
z = t.^2 ./ 2^n;
const = 2^(n/2)*pi^(n/2);

p3_true = meijerGnum(n,z) / const;
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-10 -5 0 5 10]; ax.YTick = [0 1 2 3 4];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%figure; plot(t,p3_true); title('meijerG product of 3')

figure; plot(t,p3_true - p3,'k'); %title('Error: product of 3')
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-10 -5 0 5 10]; ax.YTick = [-1.5E-15 0 1.5E-15 2.5E-15];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;


%% PLOT COEFFICIENTS of 3 GAUSSIANS

clear all; clc; %close all

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

load([path '/Runs/W010101_14dig']); % product of 3 Gaussians: W_kj
k_start = -100;
k_end = 100;
j_start = -10;
j_end = 110;
j_end = 140; % to match visual with 4 Gaussians
W_kj(122:151,:) = zeros(30,201); % to match visual with 4 Gaussians
    
figure
%surf(GMR); title('GMR matrix of coefficients')
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% surf(X,Y,log10(W_kj+eps)); %title('GMR matrix of coefficients (log scale)')
% %xlabel('translations')
% %ylabel('scales')
% %zlim([10^-10 10^-4])
% view([0 0 1])
% shading interp
% lighting phong

%USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj+eps),-16:.25:-2); %title('GMR matrix of coefficients (log scale)')
set(h,'lineStyle','none');

c=colorbar;
c.FontSize = 18;
c.Ticks=[-15:2:-3]

%title('product of 3 Gaussians')
axis([-100 100 -10 130])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-100 -50 0 50 100]; ax.YTick = [0 20 40 60 80 100 120];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
    
    
%% PLOT COEFFICIENTS of 2 GAUSSIANS

clear all; clc; %close all

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

load([path '/Runs/W0101_14dig']); % product of 2 Gaussians: W_kj
k_start = -100;
k_end = 100;
j_start = -10;
j_end = 110;
j_end = 140; % to match visual with 4 Gaussians
W_kj(122:151,:) = zeros(30,201); % to match visual with 4 Gaussians
    
figure
%surf(GMR); title('GMR matrix of coefficients')
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% surf(X,Y,log10(W_kj+eps)); %title('GMR matrix of coefficients (log scale)')
% %xlabel('translations')
% %ylabel('scales')
% %zlim([10^-10 10^-4])
% view([0 0 1])
% shading interp
% lighting phong

%USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,log10(W_kj+eps),-16:.25:-2); %title('GMR matrix of coefficients (log scale)')
set(h,'lineStyle','none');


c=colorbar;
c.FontSize = 18;
c.Ticks=[-15:2:-3];

%title('product of 2 Gaussians')
axis([-100 100 -10 130])
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-100 -50 0 50 100]; ax.YTick = [0 20 40 60 80 100 120];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% PRODUCT OF 2 GAUSSIANS: compare with besselK function
clear all; clc; close all

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
addpath( [path '/Functions'] );

load([path '/Runs/W0101_vec_14dig']); % product of 3 Gaussians: W_vec

%p3 = @(x) eval_GMR(W_vec,x);

t = linspace(-10,10,1000);
p2 = eval_GMR(W_vec,t);

figure; plot(t,p2,'k'); %title('our product of 2 Gaussians')
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-10 -5 0 5 10]; ax.YTick = [0 0.5 1 1.5];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

% %meijerG (assuming sigma_x = sigma_y = 1) same as besselK
% n = 2;
% z = t.^2 ./ 2^n;
% const = 2^(n/2)*pi^(n/2);
% 
% p2_trueG = meijerGnum(n,z) / const;
% 
% figure; plot(t,p2_trueG); title('MeijerG product of 2')
% 
% figure; plot(t,p2_trueG - p2); title('Error with MeijerG')

%besselK
s_x = 1; s_y = 1;
p2_trueK = @(x) besselk( 0, abs(x)/(s_x * s_y) ) / (pi * s_x * s_y);

figure
plot(t,p2-p2_trueK(t),'k'); %title('Error with BesselK') %error
%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-10 -5 0 5 10]; ax.YTick = [-1.5E-15 0 1.5E-15 2.5E-15];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%SET AXIS LABELS
% %labels
% ax.XTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'};
% ax.YTickLabel = {'min = -1','-0.5','0','0.5','max = 1'};


%% PLOT TWO GAUSSIANS

% f(x) is Gaussian
mu_x = 6;
s_x = 1;
f = @(x) ( 1 / (s_x*sqrt(2*pi)) ) * exp( -(x - mu_x).^2 / (2*s_x^2) );

% g(x) is Gaussian
mu_y = 2;
s_y = 1;
g = @(x) ( 1 / (s_y*sqrt(2*pi)) ) * exp( -(x - mu_y).^2 / (2*s_y^2) );

t = linspace(-20,40,10000);

figure
plot(t,f(t),'k',t,g(t),'b')

ax = gca; ax.XTick = [-20 -10 0 10 20 30 40]; ax.YTick = [0 0.1 0.2 0.3 0.4];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

    