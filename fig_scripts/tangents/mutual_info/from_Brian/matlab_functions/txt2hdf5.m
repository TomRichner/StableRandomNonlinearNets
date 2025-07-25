%%
% Convert SMR text to hdf5 format
FileName = {
    '081027a1.txt'
    '081027a2.txt'
    };
h5fn = FileName{1}(1:7);

h5name = '081027aJuxta.h5';
dsname = {
    'a1unif2'
    'a2unif2'
    };

hdf5write(h5name,dsname{1}, [],'WriteMode','overwrite');
clear data tmp
for i = 1:length(FileName)
    fid = fopen(FileName{i});
    Hdr = textscan(fid, '%s', 27);
    fn = FileName{i}(1:end-4);
    tmp = textscan(fid,'%n');
    hdf5write(h5name,dsname{i},tmp{1},'WriteMode','append');
    fclose(fid);
end
