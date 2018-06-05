%normalize
clear all; clc

load('Rns/Rn6') %load mfcc's of songs
[N,d] = size(Rn);

for i = 1:d
    mu = mean(Rn(:,i));
    st = std(Rn(:,i));
    Rn(:,1) = (Rn(:,1) - mu)/st;
end

for i=1:N
    Rn(i,:) = Rn(i,:)/norm(Rn(i,:));
end



save('Rns/Rn6n','Rn')

    