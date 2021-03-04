%Step3: Lossless comp (“gzip")
%Transmit or store with c, s & bits
M=(y.*s+2^(bits/2));
File='Compressed_B';
FMT='uint8';
fid=fopen(File,'w','l');
count=fwrite(fid,M,FMT);
status=fclose(fid);
gzip(File);