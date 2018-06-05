%project3 code - scratch sheet

clear all

load('measurement_matrix.mat')
A = measurement_matrix;
load('W2.mat')
load('W2i.mat')

img = imread('tinyhorse.tif'); 
[M, N] = size(img); 

figure
imagesc(img); colormap(gray);


% Convert image to vector
imgvec = reshape(img,M*N,1);
nlen = length(imgvec);
mlen = 400;
sigma = 0.1;


measurements = A*double(imgvec); + sigma*randn(mlen,1);

y = measurements;

epsil =50;

cvx_begin
    variable c(nlen)
    minimize( norm(c, 1));
    subject to
        norm(A*Wi*c - y,2) <= epsil;
    
cvx_end

x = Wi*c;

ximg = reshape(x,32,32);
figure
imagesc(ximg); colormap(gray);


norm(double(imgvec)-x,2)

%% TV

clear all

load('measurement_matrix.mat')
A = measurement_matrix;

img = imread('tinyhorse.tif'); 
[M, N] = size(img); 

figure
imagesc(img); colormap(gray);


% Convert image to vector
imgvec = reshape(img,M*N,1);
nlen = length(imgvec);
mlen = 400;
sigma = 0.1;


measurements = A*double(imgvec); + sigma*randn(mlen,1);

y = measurements;

epsil =1;

cvx_begin
    variable x(nlen)
    minimize( total_variation(x) );
    subject to
        norm(A*x - y,2) <= epsil;
    
cvx_end

ximg = reshape(x,32,32);
figure
imagesc(ximg); colormap(gray);

mse=sum((x-double(imgvec)).^2)/length(x)


%% TV -Rice

clear all

load('cscameradata/Phi_32.mat')

load('cscameradata/Dice_32.mat')

xt = Phi\y;

xtimg = reshape(xt,32,32);
imagesc(xtimg); colormap(gray)


nlen = 1024;
mlen = 400;
sigma = 0.1;


y = y(1:400,1);
A=Phi(1:400,:);

epsil =.0001;

cvx_begin
    variable x(nlen)
    minimize( total_variation(x) );
    subject to
        norm(A*x - y,2) <= epsil;
    
cvx_end

ximg = reshape(x,32,32);
figure
imagesc(ximg); colormap(gray);


mse=sum((x-xt).^2)/length(x)





%% check function mult_by_W


img = imread('tinyhorse.tif'); 
[M, N] = size(img);

figure
imagesc(img); colormap(gray)

% Convert image to vector
imgvec = double(reshape(img,M*N,1));

h=daubcqf(2);
level =5;

Wx = mult_by_W(imgvec,0,h,level);

Wximg = reshape(Wx,32,32);

figure
imagesc(Wximg); colormap(gray)

Wc = mult_by_W(Wx,1,h,level);

Wcimg = reshape(Wc,32,32);

figure
imagesc(Wcimg); colormap(gray)


%% create matrix W
clear all

h=daubcqf(8);
level =5;

W = zeros(1024,1024);

for k = 1:1024
    
    e = zeros(1024,1);
    e(k) = 1;
    
    wcol = mult_by_W(e,0,h,level);
    
    W(:,k) = wcol;
    
end

save('W3','W')
Wi = inv(W);
save('W3i','Wi')


%% test matrix W
clear all
load('W.mat')
load('Wi.mat')

img = imread('tinyhorse.tif'); 
[M, N] = size(img);

figure
imagesc(img); colormap(gray)

% Convert image to vector
imgvec = double(reshape(img,M*N,1));

c = W*imgvec;

cimg = reshape(c,32,32);

figure
imagesc(cimg); colormap(gray);



x = Wi*c;

ximg = reshape(x,32,32)

figure
imagesc(ximg); colormap(gray);






    
    
    
    




