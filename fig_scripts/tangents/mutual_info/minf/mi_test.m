
disp('Test: Mutual information between two images')
load mri 
A=D(:,:,8); 
B=D(:,:,9);
mi_fcn(A,A),mi_fcn(A,B)

n_shift = 100;
mi_vec = zeros(n_shift,1);
for i_shift = 1:n_shift

    B_shift=circshift(B,i_shift-1,1);

    mi_vec(i_shift,1) = mi_fcn(A,B_shift);

end

figure;
plot(mi_vec)

% disp('Test: Mutual information between two signals')
% load garchdata
% nasdaq = price2ret(NASDAQ);
% nyse = price2ret(NYSE);
% mi(nasdaq,nasdaq), mi(nasdaq,nyse)
