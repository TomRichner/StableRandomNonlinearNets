function data = choose_seq(M, N, L, rep)
%function data = choose_seq(M, N, L)
%
% M = data
% N = number of repetitions
% L = length
% rep - 0 = random without replacement; 1 = with replacement
%
% This function chooses N sequences of length L from a matrix or cell array
% of matrices M. The beginning of each sequences is choosen randomly. Each
% sequence is outputted as a column.
%
%-------------------------------------------------------------------------
% BNL, 0604

if ~iscell(M)
    tmp = M;
    clear M;
    M{1} = tmp;
end

% create data matrix, where each sequence is a column
for i = 1:length(M)
    if L > length(M{i})/4
        LL = floor(length(M{i})/4);
        disp(['Length is too long for available data in matrix ' num2str(i) '.'])
        disp(['Using L = ' num2str(LL) ' instead.'])
    else
        LL = L;
    end
    if rep == 1
        % construct table of indices for a given matrix
        tmp = floor((length(M{i})-LL-1)*rand(1,N))+1;
        index = ones(LL,1)*tmp + [0:(LL-1)]'*ones(1,N);
        data{i} = M{i}(index);
    end
    if rep == 0
        tmp = randperm(length(M{i})-LL-1);
        NN = N;
        if length(tmp)<N
            disp('CHOOSE_SEQ.M: Length of input data vector < desired repeats.')
            disp(['N = ' num2str(length(tmp)) ' for data cell ' num2str(i)]);
            NN = length(tmp);
        end
        tmp = tmp(1:NN);
        index = ones(LL,1)*tmp + [0:(LL-1)]'*ones(1,NN);
        data{i} = M{i}(index);
    end
end
