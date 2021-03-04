%% Yevgen Solodkyy
%% AMTH 308 HW4

% Uncompress
clear all
clc
close all

gunzip('BQ.gz')
file1 = 'BQ';
fid=fopen(file1,'r','l')
BQ2 = fread(fid);
status = fclose(fid);


%%
gunzip('SGN.gz')
file2 = 'SGN';
fid=fopen(file2,'r','l');
SGN2 = fread(fid);
status = fclose(fid);

SGN3=SGN2-1; % get the signs back


gunzip('Codebook.gz')
file3 = 'Codebook';
fid=fopen(file3,'r','l')
Codebook2 = fread(fid,'uint32');
status = fclose(fid);

BQ3 = Codebook2(BQ2(:)+1); % it is the codebook that is broken



BQ3=BQ3.*SGN3; 

%rearrange into a matrix
BQ3=reshape(BQ3,256,256);



%%

len = size(BQ3,1);
%build Haar filter matrix
T=[1 1; 1 -1];
I = eye(len);
H = kron(I(1:len/2,1:len/2),T)/sqrt(2); %build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);
%decode encoded image 
A=BQ3;
len2=1;
for j = 1:log2(len)
len2 = 2*len2;
P = [PT(1:len2/2,1:len2);
PB(1:len2/2,1:len2)]; 
HE = H(1:len2,1:len2);
A(1:len2,1:len2)= HE'*P'*A(1:len2,1:len2)*P*HE;
end

figure(2);
imagesc(A);
colormap gray(256);
title('Recreated')


