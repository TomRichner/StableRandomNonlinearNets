% INT CODE
% mutual information calculations
% 09/05
%
% based on np050911.m

% Find mutual information Im(del|sigma) that a single interval gives about
% the stimulus variance

clear
str = 'hst23303'; %'rdb15901'
CycleTime = 24; %sec
BinWidth = 0.02; %sec
TimeStep = 0.5; %sec time step for x-axis of plot

SpikeTimes = read_bin(str);
P = CycleTime/6;
split = [0 2*P 3*P 5*P 6*P];
% remove partial trial at end of experiment if present
Trials = round(max(SpikeTimes)/CycleTime); %complete trials
SpikeTimes = SpikeTimes(1:max(find(SpikeTimes < Trials*CycleTime)));
spikes = rem(SpikeTimes,CycleTime);

% make steady-state distributions for each of the 4 epochs
tt = 1; %width of SS epoch
for i = 1:length(split)-1
    z = split(i+1);
    tmp = diff(spikes(find(spikes > z-tt & spikes < z)));
    % get rid of negative numbers at boundaries
    tmp2 = logical((sign(tmp)+1) / 2);
    isi_ss{i} = tmp(tmp2);
end
%isi_ss{5} = [isi_ss{1};isi_ss{2};isi_ss{3};isi_ss{4}];

% combine ss distributions from epochs 2 and 4
for i = [2 4]
    z = split(i+1);
    tmp = diff(spikes(find(spikes > z-tt/2 & spikes < z)));
    % get rid of negative numbers at boundaries
    tmp2 = logical((sign(tmp)+1) / 2);
    isi_ss{i} = tmp(tmp2);
end
isi_ss{2} = [isi_ss{2}; isi_ss{4}];
% make {4} into a distribution that contains all SS ISIs
isi_ss{4} = [isi_ss{1};isi_ss{2};isi_ss{3}];

% calulate P(sigma)
for i = 1:length(isi_ss)-1
    PxSigma(i) = length(isi_ss{i})/length(isi_ss{length(isi_ss)});
end

% create probability mass functions
%mn = 0.001;
%step = 0.005;
%mx = 1;
mn = -3;
step = 0.05;
mx = 0;
dist = 'log10';
[px nx bins] = isi_to_px(isi_ss,mn,step,mx,dist,1);

% bin number for each isi
%PxBin = floor((log10(isi{5})-bins(2))./(diff(bins(2:3)))+2);
%PxBin(find(PxBin>length(px)))=length(px);


%-------------------------------------------------------------------------
% MI Calculations
% INPUTS: px nx PxSigma
%-------------------------------------------------------------------------

% MI using steady-state P(delta|sigma)
clear tmp tmp2 tmp3
sz = size(px,2);
for i = 1:size(px,2)-1
    tmp(:,i) = PxSigma(i)*px(:,i).*log2( px(:,i) ./ px(:,sz) );
    tmp2 = px(:,sz) .* log2(px(:,sz));
    tmp3(:,i) = PxSigma(i)*px(:,i).*log2( px(:,i) );
end

% MI according to: IM = sum(sum( P(s)P(r|s)log2( P(r|s)/P(r) )
mi1 = sum(sum(tmp))

% MI according to: IM = -sum(P(r)log2 P(r) + sum(sum(P(s)P(r|s)log2 P(r|s)
%mi2 = -sum(tmp2) + sum(sum(tmp3))

% where the link between these 2 equations is: P(r) = sum( P(s)P(r|s) )
% to check this:
%clear tmp4 
%tmp4 = px(:,1:3)*PxSigma';
%clf; 
%plot(tmp4); hold on
%plot(px(:,4),'r')
%legend('\Sigma P(\sigma) P(\Delta|\sigma)','P(\Delta)')


% calculate MI after a switch between sigma_est and sigma

% find maximum P for each delta
% here: px2 = P(del|sigma) and tmp = P(sig_est|del)
clear tmp tmp2
px2 = px(:,1:end-1);
nx2 = nx(:,1:end-1);

[Y I] = max(px2,[],2);
tmp = zeros(size(px2,1),size(px2,2));
for i = 1:size(px2,1)
    tmp(i,I(i)) = 1;
end
% tmp' = Psigest_delta
% px2 = Pdelta_sig

Psigest_sig = tmp'*px2;
Psigest = Psigest_sig * PxSigma';

for i = 1:size(Psigest_sig,2)
    tmp2(:,i) = PxSigma(i)*Psigest_sig(:,i).*log2(Psigest_sig(:,i)/Psigest(i));
end
mi3 = sum(sum(tmp2))


%-------------------------------------------------------------------------
% MI calculation for intervals after a switch
%-------------------------------------------------------------------------
clear tmp
for i = 1:Trials
    tmp{i} = SpikeTimes(find(SpikeTimes > (i-1)*CycleTime & SpikeTimes < i*CycleTime));
    tmp{i} = rem(tmp{i},24);
    for j = 1:length(split)-1
        a = split(j);
        z = split(j+1);
        parsed_isi{i}{j} = diff(tmp{i}(find(tmp{i} > a & tmp{i} < z)));
    end
end

% create array of first n intervals,
clear tmp
for n = 1:100
    for j = 1:length(parsed_isi{i})
        for  i = 1:length(parsed_isi)
            if n > length(parsed_isi{i}{j})
                tmp(i) = NaN;
            else
                tmp(i) = parsed_isi{i}{j}(n);
            end
        end
        tmp = tmp(~isnan(tmp));
        single_isi{n}{j} = tmp;
    end
end

% combine columns 2 and 4...
clear tmp tmp2
for i = 1:length(single_isi)
    for j = 1:size(single_isi{i},2)
        tmp{i}{j} = single_isi{i}{j};
    end
    tmp2{i}{1} = tmp{i}{1};
    tmp2{i}{2} = [tmp{i}{2} tmp{i}{4}];
    tmp2{i}{3} = tmp{i}{3};
end
single_isi = tmp2;

%get rid of any negative numbers when there aren't enough intervals in a
%given variance condition
for i = 1:length(single_isi)
    for j = 1:length(single_isi{i})
        single_isi{i}{j} = single_isi{i}{j}(find(single_isi{i}{j}>0));
    end
end

for i = 1:length(single_isi)
    [sx{i} nsx{i} bins] = isi_to_px(single_isi{i},mn,step,mx,dist,1);
end

% calculate MI after a switch between sigma_est and sigma

% find maximum P for each delta
% here: px2 = P(del|sigma) and tmp = P(sig_est|del)
for m = 1:length(sx)
    clear tmp tmp2
    [Y I] = max(sx{m},[],2);
    tmp = zeros(size(sx{m},1),size(sx{m},2));
    for i = 1:size(sx{m},1)
        tmp(i,I(i)) = 1;
    end
    %tmp' = Psigest_delta
    %sx{m} = Pdelta_sig
    Psigest_sig = tmp'*sx{m};
    
    Pdel = sum(nsx{m},2)/sum(sum(nsx{m}));
    PxSigma = sum(nsx{m}) / sum(sum(nsx{m}));
    Psigest = Psigest_sig * PxSigma';
    
    clear tmp3
    % MI (sigest | sig)
    for i = 1:size(Psigest_sig,2)
        tmp3(:,i) = PxSigma(i)*Psigest_sig(:,i).*log2(Psigest_sig(:,i)/Psigest(i));
    end
    tmp3(isnan(tmp3))=0;
    mi3(m) = sum(sum(tmp3));
    
    clear tmp4
    % MI (del | sig)
    for i = 1:size(sx{m},2)
        tmp4(:,i) = PxSigma(i)*sx{m}(:,i).*log2( sx{m}(:,i) ./ Pdel );
    end
    tmp4(isnan(tmp4))=0;
    mi4(m) = sum(sum(tmp4));
end

figure
set(gca,'FontSize',12)
plot(mi3,'.-')
xlabel('Spike interval after a variance switch')
ylabel('I_m (\sigma_{est}; \sigma) bits/interval)')
title({'Mutual Information using 3 \sigma s;'; ...
    ' ~0.16 bits/interval when using SS P(\Delta|\sigma) 1-sec distributions' ; ...
    ' ~0.41 bits/interval when using SS P(\sigma_{est}|\sigma)' })

figure
set(gca,'FontSize',12)
plot(mi4,'.-')
xlabel('Spike interval after a variance switch')
ylabel('I_m (\Delta; \sigma) bits/interval)')
title({'Mutual Information using 3 \sigma s;'; ...
    ' 0.159 bits/interval when using SS 1-sec distributions'})

