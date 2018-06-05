clear all; clc

load('Rns/Rn6n') %load mfcc's of songs
[N,n] = size(Rn);

eps = .1
eta = .01

Nt =50;



kf = @(eps,eta) (8*log( Nt*sqrt(2/eta) ) ) / eps^2;

kf(eps,eta)

k=floor( kf(eps,eta) ) +1 

K = (8*log( Nt*sqrt(2/eta) ) ) / eps^2;
K = floor(K) +1



univ_c = sqrt(k/n); %universal constant

%set = rand_pick(1,729,Nt);
%Rn = Rn(set,:);

Rn = Rn(1:Nt,:);

bad_proj_count = 0;
bp = 0;

for t = 1:5

bad_pt =[];

t;

p_dist_before = pdist(Rn)*univ_c;
 
Rk = dim_reduction(Rn,k,'s'); %reduce dimension
 
 
p_dist_after = pdist(Rk);

rel_error = p_dist_after ./ p_dist_before;

bad_pt = find(rel_error < (1-eps) | rel_error > 1/(1-eps) );
if ~isempty(bad_pt) 
    bad_proj_count = bad_proj_count +1;
    bp = bp +1;
%if bp <= 5
%subplot(5,1,bp)
end

subplot(5,1,t)

plot(rel_error,zeros(length(rel_error)),'ro') 
hold on
x = [.85:.1:1.15];
plot(x,zeros(length(x)),'k')
axis([.85 1.15 -.6 .1])
axis off
hold on
plot([1-eps 1 1/(1-eps)], zeros(3,1), 'k*' )



text('Interpreter','latex',...
	'String','$$1-\varepsilon$$',...
	'Position',[1-eps -.1],...
	'FontSize',14, ...
    'HorizontalAlignment','center')

text('Interpreter','latex',...
	'String','$$\frac{1}{1-\varepsilon}$$',...
	'Position',[1/(1-eps) -.3],...
	'FontSize',14, ...
    'HorizontalAlignment','center')

text('Interpreter','latex',...
	'String','$$1$$',...
	'Position',[1 -.2],...
	'FontSize',14, ...
    'HorizontalAlignment','center')

%end %plot
end %count


bad_proj_count
