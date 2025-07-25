function [px nx bins BinLabels] = isi_to_px(isi,mn,del,mx,dist,sub)
%function [px nx bins BinLabels] = isi_to_px(isi,mn,del,mx,dist,sub)
%
% This function takes an array of interspike interval vectors and outputs
% a histogram that is a probability mass function. This means that the
% probabilities for each x-value sum to 1. Each bin is ensured of having a
% small nonzero number.
%
% Range is in units of log10 seconds or seconds
%
% If ISI is a matrix of ISI vectors, then the vectors should be column
% vectors of ISIs.
%
% nx is total number in each bin before dividing by the total.
%
% INPUT:
% MN/MX     -   minimum/maximum range
% DEL       -   bin width for histogram
% DIST      -   type of distribution to use, 'log10' or 'linear'
% SUB       -   how zero bins are filled, Values are 0, eps, 1, or 11;
%               11 corresponds to adding 1 to every bin; default is 1;
%
% OUTPUT:
% PX        -   probability distributions
% NX        -   non-normalized distributions
% BINS      -   min edges of bins, see HISTC for details
% BinLabels -   lists the bin number corresponding to each data input
%
%-------------------------------------------------------------------------
% BNL, 05-9

% calculate probabilities
% create bin vector
bins = mn:del:mx;
bins = [-inf bins inf];

% make isi a cell array if it isn't already
if size(isi,2) > size(isi,1)
    disp('WARNING: input to ISI_to_PX may be row vectors instead of column vectors')
end

if iscell(isi)
    isi2 = isi;
else
    for i = 1:size(isi,2)
        isi2{i} = isi(:,i);
    end
end

if strcmp(dist,'log10')
    for i = 1:length(isi2)
        isi2{i} = log10(isi2{i});
    end
    dist = 'linear';
end

if strcmp(dist,'linear')
    if sub == 11
        for i = 1:length(isi2)
            [tmp BinLabels{i}] = histc((isi2{i}),bins);
            tmp = tmp + 1; %ensure each bin is nonzero;
            tmp = tmp(1:end-1); %remove +inf bin
            if tmp(1)<2 && tmp(end)<2
                tmp(1)=0; tmp(end)=0;
            else
                warning('Bin limits were surpassed; infinity bins are not empty.')
            end
            px(:,i) = tmp / sum(tmp); % normalization
            nx(:,i) = tmp;
        end
    else
        BinLabels = cell(1,length(isi2));
        for i = 1:length(isi2)
            [tmp BinLabels{i}] = histc((isi2{i}),bins);
            tmp(tmp==0) = sub; %ensure each bin is nonzero;
            tmp = tmp(1:end-1); %remove +inf bin
            if tmp(1)<2 && tmp(end)<2
                tmp(1)=0; tmp(end)=0;
            else
                warning('Bin limits were surpassed; infinity bins are not empty.')
            end
            px(:,i) = tmp / sum(tmp); % normalization
            nx(:,i) = tmp;
           
        end
    end
end

bins = bins(1:end-1); %remove +inf bin

if ~strcmp(dist,'linear')
    error('Error: distribution flag in function call is incorrect.')
end
