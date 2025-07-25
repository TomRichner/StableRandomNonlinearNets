% CONTENTS:

% PART 1: mean SNRs taken from 1000 repetitions of D for n = N
% PART 2: spikes needed to reach SNR = 1
% PART 3: time needed to reach SNR = 1
%
% BNL, 12/04




%-------------------------------------------------------------------------
% FRB217 -- Plot ISI distributions from ISI vector array
%
clear
load frb217_isi_2to10;
tmp = frb217_isi;
av = [1 30 10 6 3 10 30 60 10];
[a b] = sort(av);
isi=tmp(b);
n_isi_dist = size(isi,2); % number of ISI vectors
clear tmp

% create bin vector
mx = 2;
mn = -4;
dx = 0.05;
bins = mn:dx:mx;
bins = [-inf bins inf];

% create prob density functions
for i = 1:n_isi_dist
    tmp(i,:) = histc(log10(isi{i}),bins);
    isi_dist(i,:) = tmp(i,:)/sum(tmp(i,:));
    clear tmp;
end

% create ISI distributions for hhrk4 data
color = 'bgrcmykbgrcmyk';
mat = 1:n_isi_dist;
for i = 1:length(mat)
    plot(bins(25:120),isi_dist(mat(i),25:120),color(i))
    hold on;
end
xlabel('Log_{10} of ISI Interval')
ylabel('Probability Mass Function')
title('FRB217 ISI Distributions')
legend('Experiment 2','3','4','5','6','7','8','9','10')


%-------------------------------------------------------------------------
% HHRK4 -- Plot ISI distributions from ISI vector array
%
clear
load hhrk4_isi;
isi = hhrk4_isi;
n_isi_dist = size(isi,2); % number of ISI vectors

% create bin vector
mx = 2;
mn = -4;
dx = 0.05;
bins = mn:dx:mx;
bins = [-inf bins inf];

% create prob density functions
for i = 1:n_isi_dist
    tmp(i,:) = histc(log10(isi{i}),bins);
    isi_dist(i,:) = tmp(i,:)/sum(tmp(i,:));
    clear tmp;
end

% create ISI distributions for hhrk4 data
color = 'bgrcmykbgrcmyk';
mat = 1:n_isi_dist;
for i = 1:length(mat)
    plot(bins(25:120),isi_dist(mat(i),25:120),color(i))
    hold on;
end
xlabel('Log_{10} of ISI Interval')
ylabel('Probability Mass Function')
title('HHRK4 ISI Distributions')
legend('1','2','3','4','5','6','7','8','9','10')


%-------------------------------------------------------------------------
% FRB217 DATA
% SNR VALUE FOR A GIVEN N
%

clear
load frb217_postn1000;
aa = snr;
color = 'bgrcymk'; cl =1;
sym = '.o*';
N = [5 20 35];

for Ni = 1:3

av = [1 30 10 6 3 10 30 60 10];

% array of mean SNRs for n=N
for i = 1:size(aa,1)
    for j = 1:size(aa,2)
        pd(i,j) = aa{i,j}(N(Ni));
    end
end

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
x = round(x.*10)./10; % round to nearest tenth
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

scatter(x,yvec,[color(cl)]);
hold on;

new_x(Ni,:) = x;
new_y(Ni,:) = yvec;

% fit each side with cubic fits
xi = find(x==0);
xe = length(x);
x1 = [x(1:xi) abs(fliplr(x(1:(xi))))];
y1 = [yvec(1:xi) fliplr(yvec(1:(xi)))];
pp = polyfit(x1,y1,3);
xx = -2:.1:0;
%yy = pp(1).*xx.^5 +pp(2).*xx.^4+pp(3).*xx.^3+pp(4).*xx.^2+pp(5).*xx+pp(6);
%yy = pp(1).*xx.^4 +pp(2).*xx.^3+pp(3).*xx.^2+pp(4).*xx+pp(5);
yy = pp(1).*xx.^3+pp(2).*xx.^2+pp(3).*xx+pp(4);
plot(xx,yy,color(cl)); hold on;

x2 = [-1.*(fliplr(x(xi:xe))) x(xi:xe)];
y2 = [fliplr(yvec(xi:xe)) yvec(xi:xe)];
pp = polyfit(x2,y2,3);
xx = 0:.1:2;
%yy = pp(1).*xx.^5 +pp(2).*xx.^4+pp(3).*xx.^3+pp(4).*xx.^2+pp(5).*xx+pp(6);
%yy = pp(1).*xx.^4 +pp(2).*xx.^3+pp(3).*xx.^2+pp(4).*xx+pp(5);
yy = pp(1).*xx.^3+pp(2).*xx.^2+pp(3).*xx+pp(4);
plot(xx,yy,color(cl)); hold on;

cl = cl+1;

end

axis([-1.5 1.5 -0.5 3])
plot([-2:.1:2],ones(1,41),'k-.')
title(['FRB217 SNR vs Variance Ratios for D with n = ' int2str(N)])
xlabel(texlabel('log_{10}( sigma_{1}^2 / sigma_{2}^2)'))
ylabel('Signal-to-Noise Ratio') 
legend('5 Intervals','20 Intervals','35 Intervals')

%-------------------------------------------------------------------------
% FRB217 DATA
% Just Noticable Difference for a given spike number
%

clear
load frb217_postn1000;
aa = snr;
color = 'bgrcymk'; cl =1;
sym = '.o*';
negrat = [];
posrat = [];
warning off MATLAB:divideByZero

for N = 1:200
av = [1 30 10 6 3 10 30 60 10];

% array of mean SNRs for n=N
for i = 1:size(aa,1)
    for j = 1:size(aa,2)
        pd(i,j) = aa{i,j}(N);
    end
end

%change NaN to zeros - allows some diagonal elements to be used in analysis
pd(isnan(pd))=0; 

[r c] = find(pd); %only use nonzero elements for plotting
y = pd(find(pd));
var_row = av(r);
var_col = av(c);
avn = [1 3 6 10 30 60];

% create non-redundant matrix
for i = 1:length(avn)
    for j = 1:length(avn)
        npd(i,j) = mean(y(find(var_row==avn(i) & var_col==avn(j))));
    end
end

% use linear interpolation to find log ratio nearest SNR = 1
uptri = npd(triu(npd)>1);
uptri2 = npd(triu(npd)<1 & triu(npd)~=0);
if ~isempty(uptri)
    [tmp1 tmp2] = find(triu(npd) == min(uptri)); % value closest to 1 but > 1
    [tmp3 tmp4] = find(triu(npd) == max(uptri2)); % value closet to 1 but < 1
    yup1 = log10(avn(tmp1)./avn(tmp2));
    yup2 = log10(avn(tmp3)./avn(tmp4));
    
    m = (yup1-yup2) / (npd(tmp1,tmp2)-npd(tmp3,tmp4));
    yup = m*(1-npd(tmp3,tmp4))+yup2;
else
    yup = NaN;
end
xup = N;

dntri = npd(tril(npd)>1);
dntri2 = npd(tril(npd)<1 & tril(npd)~=0);
if ~isempty(dntri)
    [tmp1 tmp2] = find(tril(npd) == min(dntri)); % value closest to 1 but > 1
    [tmp3 tmp4] = find(tril(npd) == max(dntri2)); 
    ydn1 = log10(avn(tmp1)./avn(tmp2));
    ydn2 = log10(avn(tmp3)./avn(tmp4));
    
    m = (ydn1-ydn2) / (npd(tmp1,tmp2)-npd(tmp3,tmp4));
    ydn = m*(1-npd(tmp3,tmp4))+ydn2;
else
    ydn = NaN;
end
xdn = N;

negrat = [negrat; xup yup];
posrat = [posrat; xdn ydn];
end

% plot with exponential fits
figure;
n = 6;
nl = length(negrat);
pl = length(posrat);
x1 = negrat(4:n:nl,1);
y1 = abs(negrat(4:n:nl,2));
x2 = posrat(4:n:pl,1);
y2 = posrat(4:n:pl,2);
pp1 = polyfit(x1,log(y1),1);
yt1 = pp1(1).*x1 + pp1(2);
pp2 = polyfit(x2,log(y2),1);
yt2 = pp2(1).*x1 + pp2(2);

scatter(x1, y1, 'r'); hold on;
scatter(x2, y2,'b')
legend('Negative \sigma^{2} Ratios','Positive \sigma^{2} Ratios')
plot(x1,exp(yt1),'r')
plot(x2,exp(yt2),'b')

title('FRB217 Minimum Variance Ratio Discriminated')
xlabel('Intervals')
ylabel(texlabel('Absolute Value of Log_{10}( sigma_{1}^2 / sigma_{2}^2)'))

%-------------------------------------------------------------------------
% HHRK4 DATA
% SNR VALUE FOR A GIVEN N
%

clear
load hhrk4_postn50
aa = snr;
color = 'bgrcymk'; cl =1;
sym = '.o*';
N = [1 3 5]
for Ni = 1:3
    
av = [1 1.5871 2.5148 3.9854 6.3172 10.0028 15.8506 25.136 39.8302 63.1431];

% array of mean SNRs for n=N
for i = 1:size(aa,1)
    for j = 1:size(aa,2)
        if ~isempty(aa{i,j})
            pd(i,j) = aa{i,j}(N(Ni));
        end
    end
end
b = 4;
npd = pd([b:10],[b:10]);
[r c] = find(npd); % only use nonzero elements for plotting
x = (log10(av(r)./av(c)));
x = round(x.*10)./10; % round to nearest tenth
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

scatter(x,yvec,[color(cl)]);
hold on;
new_x(Ni,:) = x;
new_y(Ni,:) = yvec;

% fit each side with cubic fits
xi = find(x==0);
xe = length(x);
x1 = [x(1:xi) abs(fliplr(x(1:(xi))))];
y1 = [yvec(1:xi) fliplr(yvec(1:(xi)))];
pp = polyfit(x1,y1,3);
xx = -2:.1:0;
%yy = pp(1).*xx.^5 +pp(2).*xx.^4+pp(3).*xx.^3+pp(4).*xx.^2+pp(5).*xx+pp(6);
%yy = pp(1).*xx.^4 +pp(2).*xx.^3+pp(3).*xx.^2+pp(4).*xx+pp(5);
yy = pp(1).*xx.^3+pp(2).*xx.^2+pp(3).*xx+pp(4);
plot(xx,yy,color(cl)); hold on;

x2 = [-1.*(fliplr(x(xi:xe))) x(xi:xe)];
y2 = [fliplr(yvec(xi:xe)) yvec(xi:xe)];
pp = polyfit(x2,y2,3);
xx = 0:.1:2;
%yy = pp(1).*xx.^5 +pp(2).*xx.^4+pp(3).*xx.^3+pp(4).*xx.^2+pp(5).*xx+pp(6);
%yy = pp(1).*xx.^4 +pp(2).*xx.^3+pp(3).*xx.^2+pp(4).*xx+pp(5);
yy = pp(1).*xx.^3+pp(2).*xx.^2+pp(3).*xx+pp(4);
plot(xx,yy,color(cl)); hold on;

cl = cl+1;

end

axis([-1.5 1.5 -0.5 3])
plot([-2:.1:2],ones(1,41),'k-.')
title(['HHRK4 SNR vs Variance Ratios for D with n = ' int2str(N)])
xlabel(texlabel('log_{10}( sigma_{1}^2 / sigma_{2}^2)'))
ylabel('Signal-to-Noise Ratio') 
legend('1 Intervals','3 Intervals','5 Intervals')


%-------------------------------------------------------------------------
% HHRK4 DATA
% Just Noticable Difference for a given spike number
%

clear
load hhrk4_postn50
aa = snr;

negrat = [];
posrat = [];
for N = 1:20
    
av = [1 1.5871 2.5148 3.9854 6.3172 10.0028 15.8506 25.136 39.8302 63.1431];

% array of mean SNRs for n=N
for i = 1:size(aa,1)
    for j = 1:size(aa,2)
        if ~isempty(aa{i,j})
            pd(i,j) = aa{i,j}(N);
        end
    end
end
b  = 4;
npd = pd([b:10],[b:10]);
avn = av;

% use linear interpolation to find log ratio nearest SNR = 1
uptri = npd(triu(npd)>1);
uptri2 = npd(triu(npd)<1 & triu(npd)~=0);
if ~isempty(uptri)
    [tmp1 tmp2] = find(triu(npd) == min(uptri)); % value closest to 1 but > 1
    if ~isempty(uptri2)
        [tmp3 tmp4] = find(triu(npd) == max(uptri2)); % value closest to 1 but < 1
    else
        tmp3 = tmp1; tmp4 = tmp2;
    end
    yup1 = log10(avn(tmp1)./avn(tmp2));
    yup2 = log10(avn(tmp3)./avn(tmp4));
    
    m = (yup1-yup2) / (npd(tmp1,tmp2)-npd(tmp3,tmp4));
    yup = m*(1-npd(tmp3,tmp4))+yup2;
else
    yup = NaN;
end
xup = N;

dntri = npd(tril(npd)>1);
dntri2 = npd(tril(npd)<1 & tril(npd)~=0);
if ~isempty(dntri)
    [tmp1 tmp2] = find(tril(npd) == min(dntri)); % value closest to 1 but > 1
    if ~isempty(dntri2)
        [tmp3 tmp4] = find(tril(npd) == max(dntri2)); 
    else
        tmp3 = tmp1; tmp4 = tmp2;
    end
    ydn1 = log10(avn(tmp1)./avn(tmp2));
    ydn2 = log10(avn(tmp3)./avn(tmp4));
    
    m = (ydn1-ydn2) / (npd(tmp1,tmp2)-npd(tmp3,tmp4));
    ydn = m*(1-npd(tmp3,tmp4))+ydn2;
else
    ydn = NaN;
end
xdn = N;

negrat = [negrat; xup yup];
posrat = [posrat; xdn ydn];
end

% plot with exponential fits
figure;
n = 1;
nl = length(negrat);
pl = length(posrat);
x1 = negrat(1:n:nl,1);
y1 = abs(negrat(1:n:nl,2));
x2 = posrat(1:n:pl,1);
y2 = posrat(1:n:pl,2);
pp1 = polyfit(x1,log(y1),1);
yt1 = pp1(1).*x1 + pp1(2);
pp2 = polyfit(x2,log(y2),1);
yt2 = pp2(1).*x1 + pp2(2);

scatter(x1, y1, 'r'); hold on;
scatter(x2, y2,'b')
legend('Negative \sigma^{2} Ratios','Positive \sigma^{2} Ratios')
plot(x1,exp(yt1),'r')
plot(x2,exp(yt2),'b')

title('HHRK Minimum Variance Ratio Discriminated')
xlabel('Intervals')
ylabel(texlabel('Absolute Value of Log_{10}( sigma_{1}^2 / sigma_{2}^2)'))



%-------------------------------------------------------------------------
% FRB217 DATA
% SPIKES NEEDED FOR DISCRIMINATION

clear
load frb217_postn1000;
aa = snr;

av = [1 30 10 6 3 10 30 60 10];

% array of spikes needed to reach SNR = 1
for i = 1:size(aa,1)
    for j = 1:size(aa,2)
        tmp = min(find(aa{i,j}>=1));
        if isempty(tmp)
            pd(i,j) = 0;
        else
            pd(i,j) = tmp;
        end
    end
end

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

[r c] = find(npd); %only use nonzero elements for plotting
x = (log10(avn(r)./avn(c)));
x = round(x.*10)./10; % round to nearest tenth
y = npd(find(npd));

% find mean at every x-point so that only one y-value exists at each
% x-value
yvec = []; xvec = [];
for i = 1:length(x)
    yvec = [yvec mean(y(find(x==(x(i)))))];
    xvec = [xvec x(i)];
end

color = 'bgrcmykbgrcmyk';
sym = '.ox+*sd';
for i = 0:5
    plot(x(6*i+(1:6)),(y(6*i+(1:6))),[sym(i+1) '-'])
    new_mat(i+1,:) = y(6*i+(1:6));
    new_x(i+1,:) = x(6*i+(1:6));
    hold on;
end

title('Spikes Needed For Discrimination, FRB217 -- 1000 repetitions of D')
xlabel(texlabel('Log_{10}( sigma_{1}^2 / sigma_{2}^2)'))
ylabel('Number of Spikes Needed to Reach Discrimination')
legend('\sigma_{2}^2 = 1','\sigma_{2}^2 = 3','\sigma_{2}^2 = 6', ...
    '\sigma_{2}^2 = 10','\sigma_{2}^2 = 30','\sigma_{2}^2 = 60')

%-------------------------------------------------------------------------
% FRB217 DATA
% TIME NEEDED FOR DISCRIMINATION

clear
load frb217_postn1000;
load frb217_isi_2to10;
aa = snr;
bb = data;
isi = frb217_isi;

av = [1 30 10 6 3 10 30 60 10];

% array of # of spikes needed to reach SNR = 1
for i = 1:size(aa,1)
    for j = 1:size(aa,2)
        tmp = min(find(aa{i,j}>=1));
        if isempty(tmp)
            pd(i,j) = 0;
        else
            tt = cumsum(isi{i}(bb{i})); % array of times from initial seed ISI
            mt = mean(tt,2);
            pd(i,j) = mt(tmp);
        end
    end
end

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

[r c] = find(npd); %only use nonzero elements for plotting
x = (log10(avn(r)./avn(c)));
y = npd(find(npd));

% find mean at every x-point so that only one y-value exists at each
% x-value
yvec = []; xvec = [];
for i = 1:length(x)
    
    yvec = [yvec mean(y(find(x==(x(i)))))];
    xvec = [xvec x(i)];
end

color = 'bgrcmykbgrcmyk';
sym = '.ox+*sd';
for i = 0:5
    plot(x(6*i+(1:6)),(y(6*i+(1:6))),[sym(i+1) '-'])
    new_mat(i+1,:) = y(6*i+(1:6));
    new_x(i+1,:) = x(6*i+(1:6));
    hold on;
end
title('Time Needed For Discrimination, FRB217 -- 1000 repetitions of D')
xlabel(texlabel('Log_{10}( sigma_{1}^2 / sigma_{2}^2)'))
ylabel('Mean Time (sec) Needed to Reach Discrimination')
legend('\sigma_{2}^2 = 1','\sigma_{2}^2 = 3','\sigma_{2}^2 = 6', ...
    '\sigma_{2}^2 = 10','\sigma_{2}^2 = 30','\sigma_{2}^2 = 60')

%-------------------------------------------------------------------------
% HHRK4 DATA
% SPIKES NEEDED FOR DISCRIMINATION

clear
load hhrk4_postn50;
aa = snr;

av = [1 1.5871 2.5148 3.9854 6.3172 10.0028 15.8506 25.136 39.8302 63.1431];
%av = sqrt(av);

% array of spikes needed to reach SNR = 1
for i = 1:size(aa,1)
    for j = 1:size(aa,2)
        tmp = min(find(aa{i,j}>=1));
        if isempty(tmp)
            pd(i,j) = NaN;
        else
            pd(i,j) = tmp;
        end
    end
end

b = 4;
npd = pd([b:10],[b:10]);
[r c] = find(npd); %only use nonzero elements for plotting
x = (log10(av(r)./av(c)));
x = round(x.*10)./10; % round to nearest tenth
y = npd(find(npd));

% find mean at every x-point so that only one y-value exists at each
% x-value
yvec = []; xvec = [];
for i = 1:length(x)
    
    yvec = [yvec mean(y(find(x==(x(i)))))];
    xvec = [xvec x(i)];
end

color = 'bgrcmykbgrcmyk';
sym = '.ox+*sd';
for i = 0:6
    plot(x(7*i+(1:7)),y(7*i+(1:7)),[sym(i+1) '-'])
    new_mat(i+1,:) = y(7*i+(1:7));
    new_x(i+1,:) = x(7*i+(1:7));
    hold on;
end
axis([-1.3 1.3 0 38])
title('Spikes Needed For Discrimination, HHRK4 -- 1000 repetitions of D')
xlabel(texlabel('Log_{10}( sigma_{1}^2 / sigma_{2}^2)'))
ylabel('Number of Spikes Needed to Reach Discrimination')
legend('\sigma_{2}^2 = 1','\sigma_{2}^2 = 1.6','\sigma_{2}^2 = 2.5', ...
    '\sigma_{2}^2 = 4','\sigma_{2}^2 = 6.3','\sigma_{2}^2 = 10', ...
    '\sigma_{2}^2 = 15.8')



%-------------------------------------------------------------------------
% HHRK4 DATA
% TIME NEEDED FOR DISCRIMINATION

clear
load hhrk4_postn50;
load hhrk4_isi;

aa = snr;

bb = data;
isi = hhrk4_isi;

av = [1 1.5871 2.5148 3.9854 6.3172 10.0028 15.8506 25.136 39.8302 63.1431];
av = sqrt(av);

% array of # of spikes needed to reach SNR = 1
for i = 1:size(aa,1)
    for j = 1:size(aa,2)
        tmp = min(find(aa{i,j}>=1));
        if isempty(tmp)
            pd(i,j) = NaN;
        else
            tt = cumsum(isi{i}(bb{i})); % array of times from initial seed ISI
            mt = mean(tt,2);
            pd(i,j) = mt(tmp);
        end
    end
end

b=4;
npd = pd([b:10],[b:10]);
[r c] = find(npd); %only use nonzero elements for plotting
x = (log10(av(r)./av(c)));
x = round(x.*10)./10; % round to nearest tenth
y = npd(find(npd));

% find mean at every x-point so that only one y-value exists at each
% x-value
yvec = []; xvec = [];
for i = 1:length(x)
    
    yvec = [yvec mean(y(find(x==(x(i)))))];
    xvec = [xvec x(i)];
end

color = 'bgrcmykbgrcmyk';
sym = '.ox+*sd';
for i = 0:6
    plot(x(7*i+(1:7)),y(7*i+(1:7)),[sym(i+1) '-'])
    new_mat(i+1,:) = y(7*i+(1:7));
    new_x(i+1,:) = x(7*i+(1:7));
    hold on;
end

axis([-1.3 1.3 0 .8])
title('Time Needed For Discrimination, HHRK4 -- 1000 repetitions of D')
xlabel(texlabel('Log_{10}( sigma_{1}^2 / sigma_{2}^2)'))
ylabel('Mean Time (sec) Needed to Reach Discrimination')
legend('\sigma_{2}^2 = 1','\sigma_{2}^2 = 1.6','\sigma_{2}^2 = 2.5', ...
    '\sigma_{2}^2 = 4','\sigma_{2}^2 = 6.3','\sigma_{2}^2 = 10', ...
    '\sigma_{2}^2 = 15.8')
