function [y,s] = AdaptiveIntegrator(f,a,b,tol)
%returns integral y, and number of total subdivisions s

%f = @(x) log(x); a = 1; b = 2;

max_levels = 200;
stack = zeros(max_levels,2);
stack(1,:) = [a,b];

%tol = 10^-5;

I = 0;
j = 1;
max_j =1; % for keeping track of deepest level

for i = 1:10^6 
    
    
    if j == 0 
        %display('Adaptive Integrator: number of subdivisions')
        %i
        break
    end
    i;
    
    a = stack(j,1);
    b = stack(j,2);
    
    c = (a+b)/2;
    
    I1 = GaussLegendre10(f,a,b);
    I2 = GaussLegendre10(f,a,c);
    I3 = GaussLegendre10(f,c,b);
    
    if abs( I1 - (I2+I3) ) < tol
        I = I + I1;
        j = j-1;
        
    else
        stack(j,1) = c; %stack(j,2) = b;
        
        stack(j+1,1) = a; stack(j+1,2) = c;
        
        j = j+1;
        %record deepest level
        if j > max_j
            max_j =j;
        end
    
    end
    
    
    
end


y = I;
s = i;
%display max level if it is higher than 10
if max_j > 10
display(['Adaptive Integrator: max level = ' num2str(max_j)])
display(['Adaptive Integrator: total number of subdivisions= ' num2str(i)])
end

end
