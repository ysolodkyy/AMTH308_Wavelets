
% Yevgen Solodkyy
% AMTH 308 assignment 6
%% "COMPRESS"
%%
clear all
close all
clc

%read image
A0=imread('Squid.jpg');
A0=rgb2gray(A0);
A0=double(A0);
len=256;
B=imresize(A0,[len,len],'bicubic');

A1=B; % use for MSE calculations


%% Show the image

figure(1);
imagesc(B);
colormap gray(256);
title('256x256 grayscale of original')

%% build D4 filter matrix

[h0,h1,h2,h3] = solve('h0^2+h1^2+h2^2+h3^2=1', 'h0+h1+h2+h3=sqrt(2)', 'h0-h1+h2-h3=0', 'h1-2*h2+3*h3=0');

h0=(1+sqrt(3))/(4*sqrt(2));%h0(2);
h1=(3+sqrt(3))/(4*sqrt(2));%h1(2);
h2=(3-sqrt(3))/(4*sqrt(2));%h2(2);
h3=(1-sqrt(3))/(4*sqrt(2));%h3(2);

Q1 = [h0 h1;h3 -h2];
%Q1 = double(Q1)
Q2 = [h2 h3; h1 -h0];
%Q2 = double(Q2)


Q1

Q2


I = eye(len);
I2 = [I(:,len) I(:,1:len-1)];

H1 = kron(I(1:len/2,1:len/2),Q1);
H2 = kron(I2(1:len/2,1:len/2),Q2);

H = H1+H2;

% 



%build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);

%% Build the HAAR matrix
%{
Q=[1 1;1 -1];I = eye(len);
H = kron(I(1:len/2,1:len/2),Q)/sqrt(2);
% build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);
%}
%% Do full level encoding of the image

for j = 1:log2(len)-2 % added -2
    P=[PT(1:len/2, 1:len); PB(1:len/2,1:len)];
    H = H(1:len,1:len);
    H(len,1)=h1; 
    H(len-1,1)=h2;
    H(len,2)=-h0;
    H(len-1,2)=h3;
    B(1:len,1:len)=P*H*B(1:len,1:len)*H'*P'; %
    len = len/2;
end

figure(2);
image(B);
colormap gray(256); %colormap hot(16);
title('Encoded')


%% Quantize
   %cutoff = 'enter the cutoff value'
  cutoff =.9959

  len = size(B,1);
  X = sort(abs(B(:)));
  th = X(floor(cutoff*len^2));
  
 %% Built in quantization 
 bits = 8;

NP=2^(bits-1)-1; % number of partitions
[N1,N2]=size(B); % number of rows and columns -- size of the matrix
NX=N1*N2; % number of elements in the matrix?

a=abs(B(:));     % absolute value of matrix elements
SGN=sign(B(:));    % sign of matrix elements
MX=max(a);       % element of greatest magnitude

BQ = zeros(NP+1,1);
Codebook = zeros(NP+1,1);   % c stands for codebook
p = zeros(NP,1);     % p stands for partitions

Codebook(1) = 0;
d =(MX/th)^(1/NP);

for n=0:NP-1
    p(n+1)=th*d^n;
    Codebook(n+2)=th*d^(n+0);
end
BQ = quantiz(a,p);

 
 %[BQ,SGN,Codebook]= log_quant(B,th,8); 
  SGN = SGN+1; % preserve signs
 
  
 
%% TESTING -- unquantize
 

BQ3 = Codebook(BQ(:)+1);

SGN3=SGN-1; % get the signs back

BQ3=BQ3.*SGN3; 

%rearrange into a matrix
B=reshape(BQ3,256,256);


%%  Recover the image 
    % un-D4
len = size(B,1);

I2 = [I(:,len) I(:,1:len-1)];

H1 = kron(I(1:len/2,1:len/2),Q1);
H2 = kron(I2(1:len/2,1:len/2),Q2);

H = H1+H2;

%build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);


%{
% TESTING -- unhaar 

len = size(B,1);
%build Haar filter matrix
T=[1 1; 1 -1];
I = eye(len);
H = kron(I(1:len/2,1:len/2),T)/sqrt(2); %build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);
%decode encoded image 


%}

A=B;
len2=4;
for j = 1:log2(len)-2;
len2 = 2*len2;
P = [PT(1:len2/2,1:len2);
PB(1:len2/2,1:len2)]; 
HE = H(1:len2,1:len2);
    HE(len2,1)=h1; 
    HE(len2-1,1)=h2;
    HE(len2,2)=-h0;
    HE(len2-1,2)=h3;
A(1:len2,1:len2)= HE'*P'*A(1:len2,1:len2)*P*HE;
end

figure(3);
imagesc(A);
colormap gray(256);
title('Recreated')



%% calculate the PSNR 

Original_Image = A1; %double(Reference_Image);
Reconstructed_Image = A; %reshape(BQ,[256,256]); 

[M N] = size(Original_Image);
error = Original_Image - Reconstructed_Image;

error_vector=reshape(error, [256*256,1]);
MSE = sum(error_vector.*error_vector) / (M * N)

PSNR = 20*log10(255)-10*log10(MSE)

%% part 2 compress the image without encoding



%% Compress BQ
file = 'BQ';
FMT='uint8';
fid=fopen(file,'w', 'l');
count = fwrite(fid,BQ, FMT);
status = fclose(fid);
%compress the file
gzip(file);

% compression ratios
working_path= pwd;

Compressed_BQ=strcat(working_path,'\',file,'.gz');

Compressed_BQ_stats = dir(Compressed_BQ);

compressed_BQ_bytes = Compressed_BQ_stats.bytes;



%% compress the original image
file = 'A1';
FMT='uint8';
fid=fopen(file,'w', 'l');
count = fwrite(fid,A1, FMT);
status = fclose(fid);
%compress the file
gzip(file);

Compressed_A1=strcat(working_path,'\',file,'.gz');

Compressed_A1_stats = dir(Compressed_A1);

compressed_A1_bytes = Compressed_A1_stats.bytes;

'compression ratio: Quantized over original'
compression_ratio= compressed_BQ_bytes/compressed_A1_bytes
%}