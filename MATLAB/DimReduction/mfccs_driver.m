 addpath /Users/satkauskas/Documents/MATLAB/My_func/Ma/
 
 % creates mfccs, picks desired columns, stacks as one vec,
 % wrtites that vector as a row of Rn matrix
 % saves Rn matrix in Rns directory under filename
 
 %Also saves mfcc in mfcc_dir
 
 
 clear all 
 clc

 %mfcc_dir = 'MFCCs_40/';
 filename = 'Rn15'; %name of the file under which Rn is saved
 n_features = 20;

  %i_start = 10;
  %i_end = 1000; %(i_end should be less that 1653)

 %or

radius = 800;
 
 
 
 %set parameters for mfcc function
 p.fs = 11025;    %sampling freq of given wav (in Hz)
 p.visu = 0;      % figure generation
 p.hopsize = 128; % (unit: samples) aka overlap
 %p.fft_size = 256; % (unit: samples) 256 are about 23ms @ 11kHz
 p.num_ceps_coeffs = n_features;  % number of coefficents (features) ??
 p.use_first_coeff = 1;
 p.dB_max = 96; 
 
 % Load files containing filename structure
 
 load('Name_str')
 n_tracks = length(Name_str);
 
 
 min_length = 10^6;
 max_length = 0;
 j = 0; %index for checking zero v's
 
 for i = 1:n_tracks
     
     i % to see the progress
     
     load( [ 'Wavreads/' Name_str(i,:)] )
     
     mfcc = ma_mfcc(X,p);
     
     %save( [mfcc_dir Name_str(i,:) ], 'mfcc' ) 
     
     
     l = length( mfcc(1,:) );
     
         
     L(i) = l; %record the lenghts

     center = floor(l/2);
     
     % max and min length of mfcc matrix
     if  l < min_length
         min_length = l;
     end
     if l > max_length
         max_length = l;
     end
     
     %col_B = rand_pick(1,l,1000);
     %B = mfcc( :, col_B);
     
     B = mfcc( :, center-radius:center+radius);
     
     %B = mfcc( :, i_start:i_end );

     
     [n_B,m_B] = size(B);
     v = reshape(B, n_features*m_B, 1);
     
     Rn(i,:) = v;
     
     %check which v's are zero (i.e. if shit went down)
     
     if norm(v,Inf) == 0
         j=j+1;
         zero_norm(j) =i;
     end
     
         

 end
 
 
 save( ['Rns/' filename], 'Rn' )
 
 
 
 
 