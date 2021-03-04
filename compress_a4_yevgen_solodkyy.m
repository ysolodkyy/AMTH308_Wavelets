
% Yevgen Solodkyy
% AMTH 308 assignment 4
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

figure(1);
imagesc(B);
colormap gray(256);
title('256x256 grayscale of original')
%build Haar filter matrix
Q=[1 1;1 -1];I = eye(len);
H = kron(I(1:len/2,1:len/2),Q)/sqrt(2);
%build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);

%% Do full level encoding of the image
for j = 1:log2(len)
    P=[PT(1:len/2, 1:len); PB(1:len/2,1:len)];
    H = H(1:len,1:len);
    B(1:len,1:len)=P*H*B(1:len,1:len)*H'*P'; %
    len = len/2;
end

figure(2);
image(B);
colormap gray(256); %colormap hot(16);
title('Encoded')




%% Quantize
   %cutoff = 'enter the cutoff value'
  cutoff =.97

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
 
  
 %{
%% TESTING -- unquantize
 

BQ3 = Codebook(BQ(:)+1);

SGN3=SGN-1; % get the signs back

BQ3=BQ3.*SGN3; 

%rearrange into a matrix
B=reshape(BQ3,256,256);




%% TESTING -- unhaar 

len = size(B,1);
%build Haar filter matrix
T=[1 1; 1 -1];
I = eye(len);
H = kron(I(1:len/2,1:len/2),T)/sqrt(2); %build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);
%decode encoded image 
A=B;
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

%}

%% Zip BQ
file = 'BQ';
FMT='uint8';
fid=fopen(file,'w', 'l');
count = fwrite(fid,BQ, FMT);
status = fclose(fid);

%zip the file
gzip(file);

% compression ratios
working_path= pwd;

unZipped_file=strcat(working_path,'\',file);

Zipped_file=strcat(working_path,'\',file,'.gz');

unzipped_stats = dir(unZipped_file);

uncompressed_bytes = unzipped_stats.bytes;

zipped_stats = dir(Zipped_file);

compressed_bytes = zipped_stats.bytes;

BQ_compression_ratio = uncompressed_bytes/compressed_bytes

%% ZIP SGN
file2 = 'SGN';
FMT='uint8';
fid=fopen(file2,'w', 'l');
count = fwrite(fid,SGN, FMT);
status = fclose(fid);

%zip the file
gzip(file2);

%% compression ratios
working_path= pwd;

unZipped_file2=strcat(working_path,'\',file2);

Zipped_file2=strcat(working_path,'\',file2,'.gz');

unzipped_stats = dir(unZipped_file2);

uncompressed_bytes = unzipped_stats.bytes;

zipped_stats = dir(Zipped_file2);

compressed_bytes = zipped_stats.bytes;

SGN_compression_ratio = uncompressed_bytes/compressed_bytes


%% write the Codebook
file3 = 'Codebook';
FMT='uint32';
fid=fopen(file3,'w', 'l');
count = fwrite(fid,Codebook, FMT);
status = fclose(fid);
gzip(file3);
