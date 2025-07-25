function [Stim StimBit2] = make_stim(StimLength,StimMean,StimVar,SampleRate,StimBit1)
%function [Stim StimBit2] = make_stim(StimLength,StimMean,StimVar,SampleRate,StimBit1)
%
% INPUTS
% STIM_LENGTH = length of stimulus (sec)
% STIM_MEAN = desired mean of stimulus
% STIM_VAR = desired variance of stimulus
% SampleRate = sample rate of stimulus (Hz)
% STIM_BIT1 = allows priming of stimulus for continuous stimulus; is length
% of the filt
%
% OUTPUT
% STIM = 1xn vector of length STIM_LENGTH
% STIM_BIT2 = saves the last bit of the raw stimulus for later 'priming' 
%
%-------------------------------------------------------------------------
% Brian Lundstrom, Autumn 2004
% revised 4/05

if nargin == 4
    StimBit1 = [];
end

rand('state',sum(100*clock));
%randn('state',0);
noise_std = sqrt(StimVar);
half_width = 1.5*SampleRate/1000;

% EXPONENTIAL FILTER
tau = 1; %0.025; %time constant (msec)
steps = tau * SampleRate/1000; % time steps per tau
filt = exp(-(0:(10*steps))/(steps));

% GAUSSIAN FILTER
%gvar = 0.065; %variance (msec2); FWHM is ~2.35 SD, where FWHM 0.6=0.065 gvar
%filt = exp(-(-half_width:half_width).^2/(2*gvar*(SampleRate/1000)^2));

filt = filt/norm(filt);
lf = length(filt);

% explanation: 
% (1) stim = B
% (2) stimbit1 + stim = A B
% (3) conv(stim + stimbit1) = A B C
%   but conv2(...'valid') = B 
% (stimbit2 is last part of B)
if isempty(StimBit1)
    Stim = noise_std*randn(1,ceil(StimLength*SampleRate)+lf-1);
    StimBit2 = Stim(length(Stim)-lf+2:length(Stim));
    Stim = conv2(Stim,filt,'valid')+StimMean;
else
    Stim = noise_std*randn(1,ceil(StimLength*SampleRate));
    Stim = [StimBit1 Stim];
    if length(Stim) ~= (ceil(StimLength*SampleRate)+lf-1)
        error('StimBit1 does not correspond with filter length')
    end
    StimBit2 = Stim(length(Stim)-lf+2:length(Stim));
    Stim = conv2(Stim,filt,'valid')+StimMean;
end

%Stim = [0 Stim];
return;