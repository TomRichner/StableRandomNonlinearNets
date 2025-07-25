% for reading Igor files

fid=fopen(name,'r')
fseek(fid,126,-1)
data=fread(fid,inf,'float32=>double');
fclose(fid)
ll = length(data);
data = data(1:(ll-4));



% ibw files
fid=fopen(name,'r')
fseek(fid,126,-1)
data=fread(fid,inf,'int16=>int');
fclose(fid)
ll = length(data);
data = data(1:(ll-8));


data = data./40;
samp = 20;
t = (1/samp):(1/samp):300000;
st = findspikes(t,data');
