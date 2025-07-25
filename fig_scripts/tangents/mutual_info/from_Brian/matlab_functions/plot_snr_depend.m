function [new_x new_y] = plot_snr_depend(snr)
% plot independent vs. joint SNR for fly data

%clear
%load frb217_postn1000;

for n = 1:length(snr)
   
pd = snr{n};

av = [1 30 10 6 3 10 30 60 10];

%change NaN to zeros - allows some diagonal elements to be used in analysis
pd(isnan(pd))=0; 

[r c] = find(pd); %only use nonzero elements for plotting
y = pd(find(pd));
var_row = av(r);
var_col = av(c);
avn = [1 3 6 10 30 60];

% create non-reduntant matrix
for i = 1:length(avn)
    for j = 1:length(avn)
        npd(i,j) = mean(y(find(var_row==avn(i) & var_col==avn(j))));
    end
end

[r c] = find(npd);
x = log10(avn(r)./avn(c));
x = round(x.*4)./4; % round to nearest tenth
y = npd(find(npd));

yvec = [];
for i = 1:length(x)
    yvec = [yvec mean(y(find(x==(x(i)))))];
    yvec(isnan(yvec))=0;
end

[a1 a2] = unique(yvec); % get rid of repeated y values
[a3 a4] = sort(x(a2));
x = a3;
yvec = a1(a4);

new_x(:,n) = x;
new_y(:,n) = yvec;

end

plot(new_x,new_y,'.-')
title('SNR for 2 intervals, FRB217')
xlabel(texlabel('log_{10}( sigma_{2} / sigma_{1})'))
ylabel('Signal-to-Noise Ratio') 
legend('Joint','Indendent')
legend boxoff