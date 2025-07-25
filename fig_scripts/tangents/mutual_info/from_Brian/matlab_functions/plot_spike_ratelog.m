function [BinCenters,ResponsePerBin]=plot_spike_ratelog(SpikeTimes,Period,BinWidth,PlotFlag)
%function [BinCenters,ResponsePerBin]=plot_spike_rate(SpikeTimes,Period,BinWidth,Flag)
%
% plots average spike response for a repeated trial
%
% INPUT:
% SpikeTimes    =   spike time vector (sec)
% Period        =   period (sec)
% BinWidth      =   size of bin (sec)
% Flag          =   flag for displaying plots, if Flag = 0, no plots are displayed 
%
%------------------------------------------------
% BNL 4/05
Flag = PlotFlag;

%BinWidth = Period/Bins;
if nargin < 4
    PlotFlag = 0;
end
if BinWidth == 'log';
    BinFlag = 'log';
    BinWidth = 1;
end
    
% remove time before first spike if it is not in the first period
tmp = floor(SpikeTimes(1)/Period);
SpikeTimes = SpikeTimes - tmp*Period;

% delete any trailing incomplete cycle
EndTime = floor(ceil(max(SpikeTimes))/Period)*Period;
SpikeTimes = SpikeTimes(1:max(find(SpikeTimes<=EndTime)));
EndTime = ceil(SpikeTimes(end));

ss = rem(SpikeTimes,Period);
Bins = BinWidth/2:BinWidth:Period;
BW = BinWidth;

if BinFlag == 'log'
    %Bins = [0 cumsum(1.5.^(-7:.2:10))]; %Bin edges
    Bins = [0 cumsum(10.^(-1.7:.2:1.5))];
    ind = Bins<Period/2; %use only bins up to T
    Bins = [Bins(ind) Period/2]; %include the last bin, too
    Bins = [Bins Bins(2:end)+Period/2]; %defines the bin edges
    BW = diff(Bins); %binwidths
    Bins = Bins(1:end-1)+BW/2; %defines the BinCenters
end

[fr BinCenters] = hist(ss, Bins);
ResponsePerBin = fr.*1./BW.*Period./EndTime;

%{
%filter only between periods
gvar = .9^2; %variance (bins^2); FWHM is ~2.35 SD, where FWHM 0.6=0.065 gvar
half_width = 1; %ceil(100*SampleRate/1000);
filt = exp(-(-half_width:half_width).^2/(2*gvar^2));
filt = filt/norm(filt);
lf = length(filt);
lf2 = (lf-1)/2;
ResponsePerBin = conv2(ResponsePerBin,filt,'same');
%}

%{
%Pieces = 0:Period/2*SampleRate:EndTime*SampleRate;
Pieces = 0:Period/2/BinWidth:Period/BinWidth;
for i = 1:length(Pieces)-1
    tmp = ResponsePerBin(Pieces(i)+1:Pieces(i+1));
    mn = nanmean(tmp(lf2+1:end-lf2));
    tmp2 = conv2(tmp,filt,'valid'); %shorter piece is given back
    mn2 = nanmean(tmp2);
    tmp2 = tmp2/mn2*mn;
    ResponsePerBin(Pieces(i)+lf2+1:Pieces(i+1)-lf2) = tmp2;
end
%}


if Flag ~=0
    set(gca,'FontSize',16)
    plot(BinCenters,ResponsePerBin)
    grid on
    title(['Average firing rate for bin size = ' num2str(BinWidth) ' sec'])
    xlabel('Time (sec)')
    ylabel('Firing rate (Hz)')
end