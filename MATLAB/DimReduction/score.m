function [y, tie] = score(idx, D, option)

tie = 0;




%option = '1';

c1=1; c2=320; %clasical
e1=321; e2=434; % electronic
j1=435; j2=460; % jazz
m1=461; m2=505;% metal
r1=506; r2=607; % rock
w1=608; w2=729;% world



%determine each idx's genre
c_idx = find(c1<=idx & idx<=c2);
e_idx = find(e1<=idx & idx<=e2);
j_idx = find(j1<=idx & idx<=j2);
m_idx = find(m1<=idx & idx<=m2);
r_idx = find(r1<=idx & idx<=r2);
w_idx = find(w1<=idx & idx<=w2);

%determine how many of each genre you got
c_n = length(c_idx);
e_n = length(e_idx);
j_n = length(j_idx);
m_n = length(m_idx);
r_n = length(r_idx);
w_n = length(w_idx);

%determine sum of distances to each neighbor (per genre)
c_d = sum( D(c_idx) );
e_d = sum( D(e_idx) );
j_d = sum( D(j_idx) );
m_d = sum( D(m_idx) );
r_d = sum( D(r_idx) );
w_d = sum( D(w_idx) );

score_d = [c_d e_d j_d m_d r_d w_d];

if option == '1'
    
    score = [c_n e_n j_n m_n r_n w_n];
    
    y = find(score==max(score)); % finds max occurences
    
    %settling ties:
    
    if length(y) >=2
        tie = 1;
        y_idx = find( score_d(y) == min(score_d(y)) );
        if length(y_idx)>= 2
            display('double tie - BUSTED!')
        end
        y = y( y_idx );
    end
    
    
end

end






% tie = 0;
% 
% %option = '1';
% 
% c1=1; c2=256; %clasical
% e1=257; e2=347; % electronic
% j1=348; j2=368; % jazz
% m1=369; m2=404;% metal
% r1=405; r2=486; % rock
% w1=487; w2=584;% world
% 
% 
% 
% %determine each idx's genre
% c_idx = find(c1<=idx & idx<=c2);
% e_idx = find(e1<=idx & idx<=e2);
% j_idx = find(j1<=idx & idx<=j2);
% m_idx = find(m1<=idx & idx<=m2);
% r_idx = find(r1<=idx & idx<=r2);
% w_idx = find(w1<=idx & idx<=w2);
% 
% %determine how many of each gerne you got
% c_n = length(c_idx);
% e_n = length(e_idx);
% j_n = length(j_idx);
% m_n = length(m_idx);
% r_n = length(r_idx);
% w_n = length(w_idx);
% 
% %determine sum of distances to each neighbor (per genre)
% c_d = sum( D(c_idx) );
% e_d = sum( D(e_idx) );
% j_d = sum( D(j_idx) );
% m_d = sum( D(m_idx) );
% r_d = sum( D(r_idx) );
% w_d = sum( D(w_idx) );
% 
% score_d = [c_d e_d j_d m_d r_d w_d];
% 
% if option == '1'
%     
%     score = [c_n e_n j_n m_n r_n w_n];
%     
%     y = find(score==max(score)); % finds max occurences
%     
%     %settling ties:
%     
%     if length(y) >=2
%         tie = 1;
%         y_idx = find( score_d(y) == min(score_d(y)) );
%         if length(y_idx)>= 2
%             display('double tie - BUSTED!')
%         end
%         y = y( y_idx );
%     end
%     
%     
% end
% 
% end
% 


