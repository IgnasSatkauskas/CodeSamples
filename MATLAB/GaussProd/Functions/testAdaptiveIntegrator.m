%test AdaptiveIntegrator function
clear all; clc;

format long

integrand = @(x) log(x);

tol = 10^-14


tic

y=0;
    
y = AdaptiveIntegrator(integrand,0,1,tol)

[y,s] = AdaptiveIntegrator(integrand,0,1,tol)


toc




%% scratch
x=linspace(0.0001,1,1000);
plot(x,log(x),x,log(x).^3)
