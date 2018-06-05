clc


% import ground_truth.csv to get structures: VarName1 and VarName2
% i.e. names of files and their respective genres

% waveread the songs, rename them, and save in Wavereads/
% Also create char array of names

delete('Wavreads/*.*')

for i = 1:length(VarName1)
    i
    j = i;
    
    file_name = char( VarName1(i) );
    genre = char( VarName2(i) );
    genre = genre(2);
    
    if genre == 'c'
         new_name = ['c' num2str_pad(j) ];
     elseif genre == 'e'
         new_name = [ 'e' num2str_pad(j-320) ];
     elseif genre == 'j'
         new_name = [ 'j' num2str_pad(j-320-114) ];
     elseif genre == 'm'
         new_name = [ 'm' num2str_pad(j-320-114-26) ];
     elseif genre == 'r'
         new_name = [ 'r' num2str_pad(j-320-114-26-45) ];
     elseif genre == 'w'
         new_name = [ 'w' num2str_pad(j-320-114-26-45-102) ];
         
    end
    
     Name_str(i,:) = new_name;
     
    
    
    X = wavread( ['Data/' file_name(2:end-4) 'wav' ] );
    save( ['Wavreads/' new_name ] , 'X')

    
end

    save('Name_str', 'Name_str')
    
    
    