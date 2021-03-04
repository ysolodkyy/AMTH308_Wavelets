function [y,s,c] = log_quant(B,th,bits)
%Step2: Thresholding + Log Quantization
% Inputs
% B: 2-D array
% th: threshold
% bits: number of bits for quantization

% Outputs:
% y: thresholding + quantizationof abs(B)
% s: sign of B
% c: codebook

NP=2^(bits-1)-1; % number of partitions
[N1,N2]=size(B); % number of rows and columns -- size of the matrix
NX=N1*N2; % number of elements in the matrix?

a=abs(B(:));     % absolute value of matrix elements
s=sign(B(:));    % sign of matrix elements
MX=max(a);       % element of greatest magnitude

y = zeros(NP+1,1);
c = zeros(NP+1,1);   % c stands for codebook
p = zeros(NP,1);     % p stands for partitions

c(1) = 0;
d =(MX/th)^(1/NP);

for n=0:NP-1
    p(n+1)=th*d^n;
    c(n+2)=th*d^(n+0);
end
y = quantiz(a,p);