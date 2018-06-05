%import Greg's cmat dat file
V = E02;
shifts = VarName1;
scales = VarName2;

j_start = -15; j_end = 90;
k_start = -300; k_end = 300;

W = zeros(106, 601);

for n = 1:length(V);
    if ~isnan(V(n))
        j_i = scales(n) + 15 + 1;
        k_i = shifts(n) + 300 +1;
        
        W(j_i,k_i) = V(n);
        
    end
    
end

%%
clear all; clc;
load('W.mat') % load coefs (previously imported)

alpha = 0.25;
j_start = -15; j_end = 90;
k_start = -300; k_end = 300;

%Plot coefs
%% COEFFICIENTS

figure;
[X,Y]=meshgrid( k_start:k_end,j_start:j_end);

% USE COUNTOURF INSTEAD OF SURF
[~,h]=contourf(X,Y,W,-18:.25:-2); %title('GMR matrix of coefficients (log scale)')
%[~,h]=contourf(X,Y,log10(W_kj(:,:,16)+1e-32),-32:.25:-2); %title('GMR matrix of coefficients (log scale)')
set(h,'lineStyle','none');

c=colorbar;
c.FontSize = 18;
c.Ticks=[-16:2:-2];

colormap gray

%title(['Coefficients for  ' savename])
axis([k_start k_end j_start j_end])

%setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-300 -200 -100 0 100 200 300]; ax.YTick = [0 30 60 90];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;


%% plot profuct distribution
W_kj=10.^W; % un-log coefs

%% PRODUCT
tt = linspace(-20,60,5000);


p2 = eval_GMR2(W_kj,j_start,j_end,k_start,k_end,tt,alpha);
%p2 = eval_GMR2(W_kj(:,:,1),j_start,j_end,k_start,k_end,tt,alpha);

figure; plot(tt,p2,'k');

axis([-20 60 0 .08])
% %setting axis tick locations, labels and font size
ax = gca; ax.XTick = [-20 0 20 40 60]; ax.YTick = [0 .02 .04 .06 .08];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;

%% Plot factors

% f(x) is Extreme Value Distribution
mu_f = 2; % center
b_f = 3; % shape
%f = @(x) (1/b) * exp(-(x-mu)/b) .* exp( -exp(-(x-mu)/b) );
f = @(x) (1/b_f) * exp(-(x-mu_f)/b_f -exp(-(x-mu_f)/b_f));

% g(x) is a Laplace distribution
mu_g = 3;
b_g = 1;
g = @(x) (1/(2*b_g)) * exp(-abs(x-mu_g)/b_g);


t = linspace(-10,20,10000);
figure; plot(t,f(t),'k'); 
hold on; plot(t, g(t),'k'); 

%setting axis tick locations, labels and font size
axis([-10 20 0 .6])
ax = gca; ax.XTick = [-10 0 10 20]; ax.YTick = [0 .2 .4 .6];
ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;