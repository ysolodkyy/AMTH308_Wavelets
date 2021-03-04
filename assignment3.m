%% Yevgen Solodkyy
% AMTH 308
% HW #3
close all;
clear all; 
clc; 

A = ones(256);
for n=1:256
    A(n,n) =0; 
    A(1,n)=0;
    A(n,1)=0;
    A(256,n)=0;
    A(n,256)=0;
    A((257-n),n)=0;
end

for n=64:256-64
    A(64,n)=0;
    A((256-64),n)=0;
    A(n,64)=0;
    A(n,(256-64))=0;
end


A =mat2gray(A);

figure();
imshow(A)
title("matrix A")
%% construct the haar matrix
%%
[len,len] = size(A);

h = [1 1; 1 -1]; %H_1

H = 1/sqrt(2)*kron(eye(len/2),h); % H_8

%% construct the permutation matrix
%%
%len=16;
I=eye(len);

PT = I([1:2:len],:);
PB = I([2:2:len],:);

for j = 1:log2(len)
    P=[PT(1:len/2, 1:len); PB(1:len/2,1:len)];
end

%% construct B_8, 2^8x2^8: 
%%
B = P*H*A*H'*P';

figure()

subplot(2,2,1)
imshow(mat2gray(B(1:128,1:128)))
title("A (n-1)")

subplot(2,2,2)
imshow(mat2gray(B(1:128,129:256)))
title("V_n")

subplot(2,2,3)
imshow(mat2gray(B(129:256,1:128)))
title("H_n")

subplot(2,2,4)
imshow(mat2gray(B(129:256,129:256)))
title("D_n")
%% recover A
%%
%{
A_rec = H'*P'*B*P*H;

figure()
subplot(1,2,1)
imshow(A_rec)
title("original A");

subplot(1,2,2)
imshow(A_rec)
title("recovered A");
%}

%% Quantize
%%
[BQ,SGN,Codebook]= log_quant(B,0.97,8); 

file_name='BLQ_File.mat'; % 
file_id=fopen(file_name,'w');
fwrite(file_id,BQ);
fclose(file_id);
%% 
%%

%% check the compression ratio

working_path= pwd;

gzip(file_name); 

Zipped_file=strcat(working_path,'\',file_name,'.gz');

zipped_stats = dir(Zipped_file);
compressed_bytes = zipped_stats.bytes
original_bytes = 256^2;

compression_ratio = original_bytes/compressed_bytes
%%
fid=fopen(file_name,"r","l");

BQ_reopened=fread(fid);

fclose(fid);


BL = Codebook(BQ_reopened(:)+1);

% get the signs back; 
BL=BL.*SGN; 

%rearrange into matrix
BL_matrix=reshape(BL,256,256);

%max(max(BL_matrix))
figure()
imshow(mat2gray(BL_matrix))
title("B from the codebook")

%% reconstruct A
%%

A_rec = H'*P'*BL_matrix*P*H;
figure()
imshow(mat2gray(A_rec))
title("recovered A")