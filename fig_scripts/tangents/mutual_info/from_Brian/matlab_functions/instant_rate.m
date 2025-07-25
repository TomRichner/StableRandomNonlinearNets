function [Bins irate] = instant_rate(SpikeTimes,Period,BinSize)
%function [Bins irate] = instant_rate(SpikeTimes,Period,BinSize)
%
% convert a vector of spike times (sec) into an instantaneous spike rate,
% i.e. 1/ISI, at a given sampling rate
%
% BinSize   -   sec
% Period    -   average over any given period
%
%-------------------------------------------------------------------------
% B. Lundstrom 2/07
SampleRate = 500; %limiting rate to freq < 500 Hz

ind = round(SpikeTimes*SampleRate); %indices for spikes
isi = diff(ind);
if min(isi)<0
    disp('SpikeTimes vector is not monotonically increasing!')
end
tmp = zeros(1,ind(end));
ind(ind==0)=1; %round up if there is a spike at t==1;
tmp(ind) = 1;
tmp = cumsum(tmp);

ir = [zeros(1,ind(1)) SampleRate./isi(tmp(ind(1):end-1))];

%make last points and first points NaN
fillin = Period*SampleRate - rem(length(ir),Period*SampleRate);

ir(1:ind(1)) = nan;
ir = [ir nan(1,fillin)];

if ind(1)+fillin > 0.005*length(tmp)
    disp('Missing initial and final rates are > 0.5% of total.')
end

% Downsample according to BinSize:
% take mean of the data points that are being combined
scale = round(BinSize*SampleRate);
ir = ir(1:end-rem(length(ir),scale));
if scale >1
    tmp2 = nanmean(reshape(ir,scale,round(length(ir)/scale)));
else
    tmp2 = ir;
end
tmp2 = tmp2(1:end-rem(length(tmp2),Period/BinSize));
irate = reshape(tmp2,Period/BinSize,length(tmp2)/(Period/BinSize));
irate = nanmean(irate,2)';
Bins = BinSize:BinSize:Period;