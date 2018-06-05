% test function GaussLegendre10
clear all; clc;
format long

integrand = @(x) log(x);

ee = exp(1);

y = GaussLegendre10(integrand,0,1)


% for testing AdaptiveIntegrator  
% break in half the interval ...


y2 = GaussLegendre10(integrand,0,0.5) + GaussLegendre10(integrand,0.5,1)


% another half

y3 = GaussLegendre10(integrand,0,0.25) + GaussLegendre10(integrand,0.25,0.5) + GaussLegendre10(integrand,0.5,1)


y4 = GaussLegendre10(integrand,0,0.125) + GaussLegendre10(integrand,0.125,0.25) + GaussLegendre10(integrand,0.25,0.5) + GaussLegendre10(integrand,0.5,1)

y5 = GaussLegendre10(integrand,0,1/16) + GaussLegendre10(integrand,1/16,1/8) + GaussLegendre10(integrand,1/8,1/4) + GaussLegendre10(integrand,0.25,0.5) + GaussLegendre10(integrand,0.5,1)


