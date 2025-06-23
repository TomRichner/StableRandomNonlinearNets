function B = squeeze_dim(A, dim)
% SQUEEZE_DIM Remove singleton dimension from array along specified dimension
%
% B = SQUEEZE_DIM(A, DIM) returns an array B with the same elements as A 
% but with the singleton dimension DIM removed. If DIM is not a singleton 
% dimension (size > 1), the array is returned unchanged.
%
% Examples:
%   A = rand(3, 1, 4, 1, 5);
%   B = squeeze_dim(A, 2);    % Remove dimension 2, result is 3x4x1x5
%   C = squeeze_dim(A, 4);    % Remove dimension 4, result is 3x1x4x5
%
% See also SQUEEZE, PERMUTE, RESHAPE

    % Input validation
    if nargin < 2
        error('squeeze_dim:NotEnoughInputs', 'Not enough input arguments.');
    end
    
    if ~isnumeric(dim) || ~isscalar(dim) || dim < 1 || dim ~= round(dim)
        error('squeeze_dim:InvalidDimension', 'DIM must be a positive integer.');
    end
    
    % Get the size of the input array
    siz = size(A);
    
    % Check if the specified dimension exists and is singleton
    if dim > length(siz)
        % Dimension doesn't exist (implicitly size 1), return as is
        B = A;
        return;
    end
    
    if siz(dim) ~= 1
        % Dimension is not singleton, return unchanged
        B = A;
        return;
    end
    
    % Create new size vector without the specified dimension
    new_siz = siz;
    new_siz(dim) = [];
    
    % Handle edge case where we're left with empty size vector
    if isempty(new_siz)
        new_siz = 1;
    end
    
    % Reshape the array to remove the singleton dimension
    B = reshape(A, new_siz);
end