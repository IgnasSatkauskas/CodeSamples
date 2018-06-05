function [Q1, Q2, Q3, Q4, Q5] = cross_validation


G = [64 23 5 9 20 24;
     64 23 5 9 20 24;
     64 23 5 9 20 24;
     64 23 5 9 21 25;
     64 22 6 9 21 25];

q1 = G(1,:);
q2 = G(2,:);
q3 = G(3,:);
q4 = G(4,:);
q5 = G(5,:);

I = [320 434 460 505 607 729]';

Q1 = [];
Q2 = [];
Q3 = [];
Q4 = [];
Q5 = [];



s=1;

for j=1:6
    
    x = [s:I(j)];
    
    for i = 1:5
        p = randperm( length(x), G(i,j) );
        set = x(p);
        x = setdiff(x,set);
        eval(['Q' int2str(i) '=[Q' int2str(i) ' ' 'set];'] );
    end
    
    s = I(j) +1;
    
end

        
    

%%RANDOM
% c = rand_pick(1,320,64)';
% c_c = setxor(c,1:320);
% 
% e = rand_pick(321,434,23)';
% e_c = setxor(e, 321:434);
% 
% j = rand_pick(435,460,5)';
% j_c = setxor(j,435:460);
% 
% m = rand_pick(461,505,9)';
% m_c = setxor(m,461:505);
% 
% r = rand_pick(506,607,20)';
% r_c = setxor(r,506:607);
% 
% w = rand_pick(608,729,24)';
% w_c = setxor(w,608:729);
% 
% 
% Qk = Rk([c e j m r w],:);
% 
% Tk = Rk([c_c e_c j_c m_c r_c w_c],:);

%end



