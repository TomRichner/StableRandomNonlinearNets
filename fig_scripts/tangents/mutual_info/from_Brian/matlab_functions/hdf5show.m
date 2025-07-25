function name = hdf5show(filename)
% name = showhdf5(filename)
%
% show files in root directory
data = hdf5info(filename);

num = length(data.GroupHierarchy.Datasets);
for i = 1:num
    name{i,2} = i;
    name{i,1} = data.GroupHierarchy.Datasets(i).Name;
end
