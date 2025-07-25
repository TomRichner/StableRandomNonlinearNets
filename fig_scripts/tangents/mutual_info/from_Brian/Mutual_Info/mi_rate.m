% MI using rate

clear
str = 'hst23303'; %'rdb15901'
CycleTime = 24; %sec
BinWidth = 0.01; %sec
sub = eps; %what to use to ensure bins are nonzero

SpikeTimes = read_bin(str);
P = CycleTime/6;
split = [0 2*P 3*P 5*P 6*P];
% remove partial trial at end of experiment if present
Trials = round(max(SpikeTimes)/CycleTime); %complete trials
SpikeTimes = SpikeTimes(1:max(find(SpikeTimes < Trials*CycleTime)));
spikes = rem(SpikeTimes,CycleTime);

[fr BinCenters] = hist(spikes, CycleTime/BinWidth);
% ResponsePerBin = fr.*1/BinWidth.*CycleTime./EndTime;

% Epochs are 1:400 401:600 601:1000 1001:1200 for BinWidth = 0.02;

Bins = min(fr):4:max(fr);
epoch = [];
for i = 1:length(split)-1
    a = split(i)/BinWidth + 1;
    z = split(i+1)/BinWidth;
    tmp = fr(a:z);
    [tmp BC] = hist(tmp,Bins);
    tmp(find(tmp==0)) = sub; %ensure each bin is nonzero;
    Pr_sig(:,i) = tmp/sum(tmp);
end

% combine both sigma=2 distributions
Pr_sig(:,2) = (Pr_sig(:,2) + Pr_sig(:,4))/2;
Pr_sig(:,4) = [];

% calculate P_response
[tmp BC] = hist(fr,Bins);
tmp(find(tmp==0)) = sub; %ensure each bin is nonzero;
Pr = tmp'/sum(tmp);

% Calculate Steady-state and interval informations
Psig = [0.333 0.333 0.333];


% MI
clear tmp
for i = 1:size(Pr_sig,2)
    tmp(:,i) = Psig(i)*Pr_sig(:,i).*log2( Pr_sig(:,i) ./ Pr );
end

%-------------------------------------------------------------------------
% MI steady-state RATE
%-------------------------------------------------------------------------

% MI according to: IM = sum(sum( P(s)P(r|s)log2( P(r|s)/P(r) )
mi1 = sum(sum(tmp))

% steady-state --> since the distributions are non-overlapping, the SS MI
% should approach log2(choices) = log2(3) = 1.58

%-------------------------------------------------------------------------
% MI for bins after a switch
%-------------------------------------------------------------------------

% find bins
cadre = 50;
bin_num = 6;
clear Pr_sig
for i = 1:length(split)-1;
    for j = 1:bin_num
        a = split(i)/BinWidth+(j-1)*cadre+1;
        z = a+cadre-1;
        tmp = hist(fr(a:z),Bins);
        tmp(find(tmp==0)) = sub; %ensure each bin is nonzero;
        Pr_sig{j}(:,i) = tmp/sum(tmp);
     end
end
% combine both sigma=2 distributions
Pr = [];
for j = 1:bin_num
    Pr(:,j) = mean(Pr_sig{j},2);
    Pr_sig{j}(:,2) = (Pr_sig{j}(:,2) + Pr_sig{j}(:,4))/2;
    Pr_sig{j}(:,4) = [];
end


% MI (del | sig)
Psig = [0.25 .5 .25];

for j = 1:bin_num
    for i = 1:size(Pr_sig{j},2)
        tmp4(:,i) = Psig(i)*Pr_sig{j}(:,i).*log2( Pr_sig{j}(:,i) ./ Pr(:,j) );
    end
    tmp4(isnan(tmp4))=0;
    mi2(j) = sum(sum(tmp4));
end
        
