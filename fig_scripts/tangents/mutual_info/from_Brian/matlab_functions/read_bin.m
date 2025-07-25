function data = read_bin(str)
%function data = read_bin(str)
%
% str = string of filename; ".bin" will be appended
%
% BNL, 08/04

%exptype = 'hst';
%flynum = 233; 
%exp = 3;

str = [str '.bin'];
fid = fopen(str,'r','l');
data = fread(fid,'float64');
fclose(fid);