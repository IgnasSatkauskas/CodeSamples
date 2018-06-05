clear all; clc





%%CHECK DISTANCE NORMS
% load('Rns/Rn2.mat')
% K=10;
% Rk = dim_reduction(Rn,100,'s');
% [Tk,Qk] = cross_validation(Rk);
% 
% 
% 
% tree = createns(Tk);
% 
% [idx_c,D_c] = knnsearch(tree,Qk(3,:),'k',K,'distance','cityblock','IncludeTies',true);
% [idx_e,D_e] = knnsearch(tree,Qk(3,:),'k',K,'distance','euclidean','IncludeTies',true);
% 
% % [idx_c2,D] = knnsearch(Tk,Qk(3,:),'k',K,'distance','cityblock');
% % [idx_e2,D] = knnsearch(Tk,Qk(3,:),'k',K,'distance','euclidean');
% 
% idx_e = cell2mat(idx_e)
% idx_c = cell2mat(idx_c)
% 
% % idx_e2
% % idx_c2
% length( intersect(idx_e,idx_c) ) % same
