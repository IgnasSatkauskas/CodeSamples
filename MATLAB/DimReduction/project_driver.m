clear all; clc

load('Rns/Rn2') %load mfcc's of songs
[N,d] = size(Rn);
d_r = 100; %reduced dimension
univ_c = sqrt(d_r/d); %universal constant

Rn=sort(Rn);
%classifier = 'KNN'
classifier = 'LDA'


%create original index set that determines genre

Rn_idx_g = [ones(1,320) 2*ones(1,114) 3*ones(1,26) 4*ones(1,45) 5*ones(1,102) 6*ones(1,122)]';



K = 10; % # of nearest neighbors

del = [0:.2:5];
for del_i = 1:length(del) 
    
 del_i  

number_of_ties_settled = 0; %ties settled in KNN clasifier
C = zeros(6,6);

D1 = zeros(6,1);

%p_dist_before = pdist(Rn)*univ_c;

%Rn = dim_reduction(Rn,d_r,'s'); %reduce dimension


%p_dist_after = pdist(Rn);



%

for k = 1:10 % find C 10 times
    
    %Make query sets Qk
    [Q1,Q2,Q3,Q4,Q5] = cross_validation;
    
    % KNN search
    
    
    
    same_score = 0;
   
    
    for i = 1:5
        
        Q_idx = eval( ['Q' int2str(i) ] );
        % q = eval( ['q' int2str(i)] );
        T_idx = setdiff(1:729,Q_idx);
        
        Q = Rn(Q_idx,:); %query set
        T = Rn(T_idx,:); %training set
        T_idx_g = Rn_idx_g(T_idx); %corresp gernre index
        
      
        if classifier == 'KNN'
            tree = createns(T); %for KNN classifier
        elseif classifier == 'LDA'
            [M,s,S,Pi,drop_coor] = LDA(T,T_idx_g,del(del_i)); %for LDA clasifier
        end
        
        
        
        
       
        %CLASSIFY
        for j=1:length(Q_idx)
            
            
            if classifier == 'KNN'  % KNN Classifier
                [idx,D] = knnsearch(tree,Q(j,:),'k',K,'distance','euclidean');
                
                id = T_idx(idx);
                
                [m,tie] = score(id,D,'1'); % clasifies given song
                number_of_ties_settled = number_of_ties_settled + tie;
                
            
            elseif classifier == 'LDA' % LDA Classifier
                [m,tie] = LDA_score(M,s,S,Pi, Q(j,:) );
                
                number_of_ties_settled = number_of_ties_settled + tie;
                
            end
            
            
           
            % Confusion matrix
            id_j=Q_idx(j); %index of jth query pt
            
            n = Rn_idx_g( id_j); %true genre
            
            
            D1(n) = D1(n)+1;
            
            C(n,m) = C(n,m) +1;
            
            
            
            
            
        end %j=1:length(Q_idx)
        
        
        
        
        
    end %i=1:5
    
%     for g=1:6
%         figure(g)
%         plot(sort(D1(g,:)))
%         axis([0 10 0 300])
%     end
    
    
    
    C;
    sumC = sum(sum(C));
    
    
end %for k=1:10

C = C/10;
correct_proc_overall = sum(diag(C))/729;
diag_true = [320 114 26 45 102 122]';
correct_proc_norm = mean( diag(C)./ diag_true );

number_of_ties_settled;

c_o(del_i) = correct_proc_overall;
c_n(del_i) = correct_proc_norm;
c_d(del_i) = drop_coor;
end %for del_1 = 1: length(del)

subplot(2,1,1)
plot(del,c_o,'ko',del,c_n,'bo')
legend('overall','normalized')
xlabel('\Delta')
ylabel('Accuracy')
subplot(2,1,2)
plot(del,c_d,'ro')
ylabel('# of dropped coord')
xlabel('\Delta')

%%TRASH?
%CLASIFIER

% %score
% c = length( find(idx<=256) )/K; %clasical
% e = length( find(idx>=257 & idx<=347) )/K; % electronic
% j = length( find(idx>=348 & idx<=368) )/K; % jazz
% m = length( find(idx>=369 & idx<=404) )/K; % metal
% r = length( find(idx>=405 & idx<=486) )/K; % rock
% w = length( find(idx>=487 & idx<=584) )/K; % world
%
% score = [c e j m r w];
%
% s_inx = find( score == max(score) );
%
% %check if there is more that one max scor
% if length(s_inx) >= 2
%    same_score = same_score +1;
% end
%
% m = s_inx(1:1); %if there is more than one max score pick first
%
