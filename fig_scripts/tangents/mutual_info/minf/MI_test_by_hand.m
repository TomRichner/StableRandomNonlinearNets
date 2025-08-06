


fs = 100;
dt = 1/fs;
t = dt:dt:100;
f = 2;
X = sin(2*pi*2*t)';
Y = sin(2*pi*2*t)';

n_bins = 16;
MI_mi_fcn = mi_fcn(X,Y,n_bins+1)

edges = linspace(-1,1,n_bins+1);
p_x = histcounts(X,edges);
p_y = histcounts(Y,edges);

p_joint = zeros(n_bins,n_bins);

X_levels = (X+1)/2*

for i_x = 1:n_bins
    for i_y = 1:n_bins



        p_joint(i_x,i_y) = 


    end
end

