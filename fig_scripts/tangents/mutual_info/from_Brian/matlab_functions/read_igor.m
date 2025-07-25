function data = read_igor(filename, type)
% function data = read_igor(filename, type)
% read files from Igor binaries such as .bwav
%
% INPUT
% type      = 'int16' or 'float32'

% ibw files
%fid=fopen(filename,'r');
fid=fopen(filename,'r','l');

%skips the 126 bytes of header
fseek(fid,126,-1);
data=fread(fid,inf,type);
fclose(fid);

%take off last 16 bytes
if(strcmp(type,'int16'))
    data=data(1:length(data)-8);
elseif(strcmp(type,'float32'))
    data=data(1:length(data)-4);
end

