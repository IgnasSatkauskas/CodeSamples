%scratch

%% time score
% 
% t=linspace(0,60,300); %in min
% time_score = t.^2./100;
% plot(t,time_score)
% 
% %conclusion - under 20 min = less than 5pt

%% check wavelets scores
clear all

load('score_sym.mat')
load('img_names.mat')

for k = 1:10
    
    img = imread(['grayscale256/' img_names{k}]);
    figure(k)
    imagesc(img); colormap(gray);
    figure(k+10*k)
    surf(score(3:8,:,k))
end

%% test one image wavelets
clear all;

img_names = {'Splash.tiff','Plane.tiff','Peppers.tiff','Man.tiff',...
            'Lena.tiff','Goldhill.tiff', 'Cameraman.tiff','Boat.tiff','Barbara.tiff','Baboon.tiff'};

    
img = imread(['grayscale256/' img_names{8}]);
%img = left;
%figure; imagesc(img);  colormap(gray);

n_h = 10;
n_levels =6;
n_coef = 2000;

[coefficients locations] = test_ImageEncode(img,n_h,n_levels,n_coef); %#ok<*NCOMMA>
            
recon_img = test_ImageDecode(coefficients, locations,n_h,n_levels,n_coef);

%sum_img = recon_img + recon_curvlet_img;

            
score = calcMSE(img,recon_img)

figure; imagesc(recon_img);  colormap(gray);
figure; imagesc(img);  colormap(gray);



 %% Testing # of coeficients
clear all;

img_names = {'Splash.tiff','Plane.tiff','Peppers.tiff','Man.tiff',...
            'Lena.tiff','Goldhill.tiff', 'Cameraman.tiff','Boat.tiff','Barbara.tiff','Baboon.tiff'};

img = imread(['grayscale256/' img_names{7}]);
%figure; imagesc(img);  colormap(gray);

n_h = 2;
n_levels =3;
 
h = daubcqf(n_h); % wavelet must be even

w = Idwt2(double(img),h,n_levels); %wavelet coeficients 

figure; imagesc(w);  colormap(gray);

wreshaped = reshape(w, 256^2, 1); 

 %% testing curvelet package : transform
 clear all;
 
 addpath(genpath('/Users/satkauskas/Documents/MATLAB/My_tool_boxes/CurveLab-2.1.2/fdct_usfft_matlab'));
 
 img_names = {'Splash.tiff','Plane.tiff','Peppers.tiff','Man.tiff',...
            'Lena.tiff','Goldhill.tiff', 'Cameraman.tiff','Boat.tiff','Barbara.tiff','Baboon.tiff'};
        
 img = imread(['grayscale256/' img_names{7}]);
 
 n_levels = 8; % 3 scales: max is 8 (2^8=256)
 
 C = fdct_usfft(double(img),1,n_levels); %find curvlet transform
 
 %low pass - 1st element of C
 figure;
 imagesc(C{1}{1}); colormap(gray); 
 
 %high pass - last element of C
 figure;
 imagesc(C{end}{1}); colormap(gray);
 
 %diferent orientation curvelets - middle elements of C
 figure;
 
 scale = 2; %goes from 2 to 7 (from low to hight)
 n_orient = 7; % up to 8
 imagesc(C{scale}{n_orient}); colormap(gray);
 
 figure;
 scale = 7;
 n_orient = 20; % up to 64
 imagesc(C{scale}{n_orient}); colormap(gray);
 
 
 %% testing curvelet package: inverse transform
 clear all
 is_real = 1;
 
 img_names = {'Splash.tiff','Plane.tiff','Peppers.tiff','Man.tiff',...
            'Lena.tiff','Goldhill.tiff', 'Cameraman.tiff','Boat.tiff','Barbara.tiff','Baboon.tiff'};
        
 img = imread(['grayscale256/' img_names{2}]);
 
 n_levels = 2;
 
 C = fdct_usfft(double(img),is_real); %find curvlet transform
 
 %D = {C{1} C{7} C{end}};
 
 recon_img = ifdct_usfft(C,is_real); %inverse
 
 score = calcMSE(img,recon_img)
 
 
 
 figure; imagesc(recon_img); colormap(gray);
 
 
 
 
 %% testing drop smallest curvelet coef
 
 clear all
 is_real = 1;
 n_levels = 3;
 
 img_names = {'Splash.tiff','Plane.tiff','Peppers.tiff','Man.tiff',...
            'Lena.tiff','Goldhill.tiff', 'Cameraman.tiff','Boat.tiff','Barbara.tiff','Baboon.tiff'};
 
        
 
     
 img = imread(['grayscale256/' img_names{10}]);
 
 
 C = fdct_usfft(double(img),is_real,n_levels); %find curvlet transform
 
 %drop low and high pass
 
 C{1}{1}=zeros(size(C{1}{1}));
 C{end}{1}=zeros(size(C{end}{1}));
 
 
 
 %find drop value of the coeficient
 
 coef =[];
 
 for k = 1:length(C)
     for m = 1:length(C{k})
         coef = [coef; C{k}{m}(:)]; % stack coef in one vec
     end
 end
 
 %sort coefficients
 length_coef=length(coef);
 coef = abs(coef);
 coef = sort(coef,'descend');
 n_coef =1000;
 drop_value = coef(n_coef);
 
 %zero out coef below drop value
 for k = 1:length(C)
     for m = 1:length(C{k})
         C{k}{m} = C{k}{m}.* (abs(C{k}{m})>drop_value);
     end
 end
 
 
 percentage_dropped = (length_coef-n_coef)/length_coef;
 
 
 recon_img = ifdct_usfft(C,is_real); %inverse
 
 score = calcMSE(img,real(recon_img))
 
 
 
 figure; imagesc(real(recon_img)); colormap(gray);

 figure; imagesc(real(img)); colormap(gray);

 left = double(img)-recon_img;
 recon_curvlet_img = recon_img;
 
 
 %% Cell business 
clear all;
clc

A = {5, 'go', 'on a'} %1x3 cell
B = {'dick',[1 2; 3 4]} %1x2 cell
C = {A,B} % 1x2 cell

%access matrix

C{1,2}{1,2}
%or 
C{2}{2}

%write 'go on a dick'
[C{1,1}{1,2} ' ' C{1,1}{1,3} ' ' C{1,2}{1,1} ]
%or
[C{1}{2} ' ' C{1}{3} ' ' C{2}{1}]


%% Lexicographic business
clear all;
clc

al = 8
ak = 10
A=rand(ak,al);
A=abs(A);

k=10;l=8;
A(k,l)=0;
A

ind=find(A==min(min(A)))

if ind-ak*floor(ind/ak)==0 
    gk=ak
    gl=floor(ind/ak)
else
    
    gk=ind - ak*floor(ind/ak)
    gl=floor(ind/ak)+1
end


%%again lexicographic - Im super idiot!
A = rand(10,10)
k=10;
l=3;
A(k,l)=0;
A

[m1,ind1]=min(A)

[m2,ind2]=min(m1)

k
l
gk=ind1(ind2)
gl=ind2

%% test matlab and rice wavelet filters
n_hm = 1; % scaling matlab
n_hr = n_hm*2; %scaling rice


hr = daubcqf(n_hr) % already normalized

hm = dbaux(n_hm)./norm(dbaux(n_hm)) % + normalization



%% plot Db total scores
total_score = zeros(8,20);
for i=1:10
    total_score =  total_score + score(:,:,i)
end


surf(1:20,3:8,total_score(3:8,:))
ylabel('wavelet scale')
xlabel('wavelet filter')
zlabel('Total RMS error')

figure;
score_cam = score(:,:,7);

surf(1:20,3:8,score_cam(3:8,:))
ylabel('wavelet scale')
xlabel('wavelet filter')
zlabel('Total RMS error')

%axis([1 20 3 8 1100 1500])




    