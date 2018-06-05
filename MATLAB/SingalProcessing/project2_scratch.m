
%%First atempt

clear all;

%load data

load ROITimeseries
load ratings

%create training set from 1st video:combine L and R 

ROI_tr = [ROIsLHvideo1 ROIsRHvideo1];

ratings_tr = video1ratings;

%create test set from the 1st half of the second video

ROI_ts = [ROIsLHvideo2(1:434,:) ROIsRHvideo2(1:434,:)];

ratings_ts_true = video2ratingshalf;


%simple regression

% do all ratings

% A = ROI_tr;
% B = ratings_tr;

% introduce lag 

lag = 3;
A = ROI_tr(lag+1:end,:);
B = ratings_tr(1:end-lag,:);


Alpha = A\B;


%predict ratings of test set

ratings_ts = zeros(434,30);

for k =1:30
    
    ratings_ts(:,k) = ROI_ts*Alpha(:,k);
    
end

% %space business
% space = [ratings_ts(:,14:19) ratings_ts(:,30)];
% 
% for ka =1:size(space,1)
%     a = space(ka,:);
%     indx = find(a==max(a));
%     z = zeros(1,length(a));
%     z(indx) = 1;
%     a = a.*z;
%     a(indx) = 1;
%     space(ka,:) = a;
%     
%     
% end
% 
% 
% ratings_ts(:,14:19) = space(:,1:6);
% ratings_ts(:,30) = space(:,7);



%save stuf under right names
solution = ratings_ts; %#ok<NASGU>
trueratings = ratings_ts_true;

save('FMRIMindReadingSolution','solution')
save('TrueRatings','trueratings')   





%% Now make training set from 1.5 movies for submission


clear all;

%load data

load ROITimeseries
load ratings

%create training set from 1st video:combine L and R 

ROI_tr = [ROIsLHvideo1 ROIsRHvideo1];

ratings_tr = video1ratings;

%create test set from the 1st half of the second video

ROI_ts = [ROIsLHvideo2(1:434,:) ROIsRHvideo2(1:434,:)];

ratings_ts_true = video2ratingshalf;

%full training set for submission
ROI_full = [ROI_tr; ROI_ts];
ratings_full = [ratings_tr; ratings_ts_true];

%test set for submission
ROI_2ndhalf = [ROIsLHvideo2(435:end,:) ROIsRHvideo2(435:end,:)];

%simple regression

% do all ratings

% A = ROI_tr;
% B = ratings_tr;

% introduce lag 

lag = 3;
A = ROI_full(lag+1:end,:);
B = ratings_full(1:end-lag,:);


Alpha = A\B;


%predict ratings of last half

ratings_2ndhalf = zeros(434,30);

for k =1:30
    
    ratings_2ndhalf(:,k) = ROI_2ndhalf*Alpha(:,k);
    
end


%save stuf under right names

solution = ratings_2ndhalf;

save('FMRIMindReadingSolution','solution')

%% Second atempt : PCA find best dim

clear all;

%load data

load ROITimeseries
load ratings

%first reduce dimention of entire set i.e. both movies

ROI1 = [ROIsLHvideo1 ROIsRHvideo1].'; 
ROI2 = [ROIsLHvideo2 ROIsRHvideo2].'; 
X = [ROI1 ROI2];

%X=[1 2 3 4; 2 0 0 1; 1 1 1 1];

%center data

n=size(X,2);
C = eye(n) - ones(n,1)*ones(1,n)/n;
X = X*C;


% %kernel 
% ker= @(x,y) dot(x,y);
% 
% K = zeros(n,n);
% 
% for i=1:n
%     for j=1:n
%         K(i,j) = ker(X(:,i),X(:,j));
%     end
% end
% 
% [V,D] = eig(K);


K = X*X.';

[V,D]=eig(K);

d=diag(D);


d=d(end:-1:1);
V=V(:,end:-1:1);

for kd = 1:60
    kd

dim = kd*5;

U = V(:,1:dim);

Xnew = zeros(dim,n);
for k =1:n
    Xnew(:,k) = U.'*X(:,k);

end

Xnew=Xnew.';

%create training set 

ROI_tr = Xnew(1:858,:);

ratings_tr = video1ratings;

%create test set

ROI_ts = Xnew(859:1292,:);

ratings_ts_true = video2ratingshalf;


%simple regression

% do all ratings

% A = ROI_tr;
% B = ratings_tr;

% introduce lag 

lag = 3;
A = ROI_tr(lag+1:end,:);
B = ratings_tr(1:end-lag,:);


Alpha = A\B;


%predict ratings of test set

ratings_ts = zeros(434,30);

for k =1:30
    
    ratings_ts(:,k) = ROI_ts*Alpha(:,k);
    
end


%save stuf under right names
solution = ratings_ts; 
trueratings = ratings_ts_true;

save('FMRIMindReadingSolution','solution')
save('TrueRatings','trueratings')   

%SCORE

load FMRIMindReadingSolution.mat
load TrueRatings.mat

% Initialize weighting vector.  
weights = [3*ones(1,13) ones(1,17)];

%Compute correlation between true and estimated ratings for each of 30 features.  
for k = 1:30
    C = cov(trueratings(:,k),solution(:,k),1); 
    if (C(1,1) ~= 0 && C(2,2) ~= 0) % neither vector is constant
        R(k) = C(1,2)/sqrt(C(1,1)*C(2,2));
    else 
        R(k) = 0;
        
    end
end
  
 
% Score is weighted sum of correlations.  
score(kd) = sum(weights.*R); 
disp(sprintf('Student score is: %d', score))

end

figure
plot(5:5:300,score,'ro')
xlabel('Dimension d')
ylabel('Score')
%title('Score')

figure
plot(log10(d),'go')
xlabel('Dimension d')
ylabel('Eigenvalue')
%title('Eigs')

%% simple PCA
clear all;

%load data

load ROITimeseries
load ratings

%first reduce dimention of entire set i.e. both movies

ROI1 = [ROIsLHvideo1 ROIsRHvideo1].'; 
ROI2 = [ROIsLHvideo2 ROIsRHvideo2].'; 
X = [ROI1 ROI2];

%X=[1 2 3 4; 2 0 0 1; 1 1 1 1];

%center data

n=size(X,2);
C = eye(n) - ones(n,1)*ones(1,n)/n;
X = X*C;


% %kernel 
% ker= @(x,y) dot(x,y);
% 
% K = zeros(n,n);
% 
% for i=1:n
%     for j=1:n
%         K(i,j) = ker(X(:,i),X(:,j));
%     end
% end
% 
% [V,D] = eig(K);


K = X*X.';

[V,D]=eig(K);

d=diag(D);

d=d(end:-1:1);
V=V(:,end:-1:1);


dim = 17;

U = V(:,1:dim);

Xnew = zeros(dim,n);
for k =1:n
    Xnew(:,k) = U.'*X(:,k);

end

Xnew=Xnew.';

%create training set 

ROI_tr = Xnew(1:858,:);

ratings_tr = video1ratings;

%create test set

ROI_ts = Xnew(859:1292,:);

ratings_ts_true = video2ratingshalf;


%simple regression

% do all ratings

% A = ROI_tr;
% B = ratings_tr;

% introduce lag 

lag = 3;
A = ROI_tr(lag+1:end,:);
B = ratings_tr(1:end-lag,:);


Alpha = A\B;


%predict ratings of test set

ratings_ts = zeros(434,30);

for k =1:30
    
    ratings_ts(:,k) = ROI_ts*Alpha(:,k);
    
end




ratings_ts(:,14:19) = space(:,1:6);
ratings_ts(:,30) = space(:,7);



%save stuf under right names
solution = ratings_ts; 
trueratings = ratings_ts_true;

save('FMRIMindReadingSolution','solution')
save('TrueRatings','trueratings')   


%% Kernel PCA

tic

clear all;

%load data

load ROITimeseries
load ratings

%first reduce dimension of entire set i.e. both movies

ROI1 = [ROIsLHvideo1 ROIsRHvideo1].'; 
ROI2 = [ROIsLHvideo2 ROIsRHvideo2].'; 
X = [ROI1 ROI2];

n=size(X,2);

%center data

C = eye(n) - ones(n,1)*ones(1,n)/n;
X = X*C;


% %kernel 
% ker= @(x,y) dot(x,y);

%gaussian kernel
s = 700;
ker = @(x,y) exp(-dot(x-y,x-y)/(2*s^2));


K = zeros(n,n);

for i=1:n
    for j=1:n
        K(i,j) = ker(X(:,i),X(:,j));
    end
end



[V,D] = eig(K); %find eigs


d=diag(D);

d=d(end:-1:1);
V=V(:,end:-1:1);



dim = 18;

U = V(:,1:dim);

%normalize

for kd=1:dim
    u = U(:,kd);
    un = u / sqrt((u.'*K*u));
    U(:,kd)=un;
end


Xnew = U.'*K;



Xnew=Xnew.';

%create training set 

ROI_tr = Xnew(1:858,:);

ratings_tr = video1ratings;

%create test set

ROI_ts = Xnew(859:1292,:);

ratings_ts_true = video2ratingshalf;


%simple regression

% do all ratings

% A = ROI_tr;
% B = ratings_tr;

% introduce lag 

lag = 3;
A = ROI_tr(lag+1:end,:);
B = ratings_tr(1:end-lag,:);


Alpha = A\B;


%predict ratings of test set

ratings_ts = zeros(434,30);

for k =1:30
    
    ratings_ts(:,k) = ROI_ts*Alpha(:,k);
    
end


%save stuf under right names
solution = ratings_ts; 
trueratings = ratings_ts_true;

save('FMRIMindReadingSolution','solution')
save('TrueRatings','trueratings')   

%SCORE

load FMRIMindReadingSolution.mat
load TrueRatings.mat

% Initialize weighting vector.  
weights = [3*ones(1,13) ones(1,17)];

%Compute correlation between true and estimated ratings for each of 30 features.  
for k = 1:30
    C = cov(trueratings(:,k),solution(:,k),1); 
    if (C(1,1) ~= 0 && C(2,2) ~= 0) % neither vector is constant
        R(k) = C(1,2)/sqrt(C(1,1)*C(2,2));
    else 
        R(k) = 0;
        
    end
end
  
 
% Score is weighted sum of correlations.  
score = sum(weights.*R); 
disp(sprintf('Student score is: %d', score))

toc


%% Kernel PCA for submission

tic

clear all;

%load data

load ROITimeseries
load ratings

%first reduce dimention of entire set i.e. both movies

ROI1 = [ROIsLHvideo1 ROIsRHvideo1].'; 
ROI2 = [ROIsLHvideo2 ROIsRHvideo2].'; 
X = [ROI1 ROI2];

n=size(X,2);

%center data

C = eye(n) - ones(n,1)*ones(1,n)/n;
X = X*C;


% %kernel 
% ker= @(x,y) dot(x,y);

%gaussian kernel
s = 700;
ker = @(x,y) exp(-dot(x-y,x-y)/(2*s^2));


K = zeros(n,n);

for i=1:n
    for j=1:n
        K(i,j) = ker(X(:,i),X(:,j));
    end
end


[V,D] = eig(K); %find eigenvectors


d=diag(D);

%sort from largest to smallest
d=d(end:-1:1);
V=V(:,end:-1:1);


dim = 18; %pick # of eigenvectors

U = V(:,1:dim);

%normalize eigenvectors

for kd=1:dim
    u = U(:,kd);
    un = u / sqrt((u.'*K*u));
    U(:,kd)=un;
end


Xnew = U.'*K; %new points (reduced dimension)

Xnew=Xnew.';

%create training set 

ROI_tr = Xnew(1:858,:);

ratings_tr = video1ratings;

%create test set

ROI_ts = Xnew(859:1292,:);

ratings_ts_true = video2ratingshalf;


%full training set for submission
ROI_full = [ROI_tr; ROI_ts];
ratings_full = [ratings_tr; ratings_ts_true];

%test set for submission
ROI_2ndhalf = Xnew(1293:end,:);



%simple regression

% introduce lag 

lag = 3;
A = ROI_tr(lag+1:end,:);
B = ratings_tr(1:end-lag,:);


Alpha = A\B;


%predict ratings of test set

ratings_2ndhalf = zeros(434,30);

for k =1:30
    
    ratings_2ndhalf(:,k) = ROI_2ndhalf*Alpha(:,k);
    
end


%save stuf under right names
solution = ratings_2ndhalf; 

save('FMRIMindReadingSolution','solution')


% bingo!

