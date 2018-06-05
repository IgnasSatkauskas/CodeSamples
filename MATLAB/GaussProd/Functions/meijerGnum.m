function y = meijerGnum(n,t)
% INPUT: t vector of points where meijerG is evaluated
% OUTPUT: y values of meijerG_2020(0,0|t)
% OUTPUT: y values of meijerG_3030(0,0,0|t)
% OUTPUT: y values of meijerG_4040(0,0,0,0|t)

% n = 2, n = 3 or n = 4 only

%n product of n Gaussians 
% y is the same size as t
% meijerG parameters are "for product of n=2,3,4 gaussians"

% evaluate at t^2 / 2^n*(sigma_x*sigma_y*sigma_z)^2
% and devide by 2^(n/2)*pi^(n/2)*sigma_x*sigma_y*sigma_z

precision = 16;
y = zeros(size(t));

for j=1:length(t)
    
    
    z = num2str(t(j),precision);
    if (n == 2)
    string = ['meijerG([[],[]],[[0,0],[]],' z ')'];
    end
    if (n == 3)
    string = ['meijerG([[],[]],[[0,0,0],[]],' z ')'];
    end
    if (n == 4)
    string = ['meijerG([[],[]],[[0,0,0,0],[]],' z ')'];
    end
    
    %y = evalin(symengine,'meijerG([[],[]],[[0,0,0],[]],z)');
    g = evalin(symengine,string);
    y(j)=double(g);
    
end



