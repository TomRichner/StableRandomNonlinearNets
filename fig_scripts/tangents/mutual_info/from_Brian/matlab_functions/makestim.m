function [Stim StimBit2] = MakeStim(StimLength,StimMean,StimVar,SampleRate,StimBit1)
% [stim] = makestim(stim_length)
%
% INPUTS
% STIM_LENGTH = length of stimulus (sec)
% STIM_MEAN = desired mean of stimulus
% STIM_VAR = desired variance of stimulus
% SampleRate = sample rate of stimulus (Hz)
% STIM_BIT1 = allows priming of stimulus for continuous stimulus; is length
% of the filt, i.e. length 61 for the Gaussian filter.
%
% OUTPUT
% STIM = 1xn vector of length STIM_LENGTH
% STIM_BIT2 = saves the last bit of the raw stimulus for later 'priming' 
%
%-------------------------------------------------------------------------
% Brian Lundstrom, Autumn 2004
% revised 4/05

rand('state',sum(100*clock));

noise_std = sqrt(StimVar);

% EXPONENTIAL FILTER
%filt = exp(-(0:30)/10);

% GAUSSIAN FILTER
filt = exp(-(-30:30).^2/100);
filt = filt/norm(filt);
lf = length(filt);

Stim = noise_std*randn(1,StimLength*SampleRate+lf-1);
StimBit2 = Stim(length(Stim)-lf+2:length(Stim));

Stim = [StimBit1 Stim];
Stim = conv(Stim,filt)+StimMean;
if isempty(StimBit1)
    Stim = Stim(1:(length(Stim)-lf+1));
else
    Stim = Stim(lf+1:(length(Stim)-lf+1));
end

return;