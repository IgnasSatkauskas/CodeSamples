%command line test file

clear all;
p = 0;
for i=1:20
    i
    p = p + 1;
end

path = '/Users/satkauskas/Documents/MATLAB/GaussProdProject/New';
save([path '/Runs/p'],'p');
p

