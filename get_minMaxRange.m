function [min_max_range] = get_minMaxRange(n,n_a,n_b)

    % check if all inputs are integer and non negative.  n_states must be 1 or greater, n must be 1 or greater. n_a and n_b can be 0.
    
    % Check n
    if ~isnumeric(n) || ~isscalar(n) || floor(n) ~= n || n < 1
        error('n must be a positive integer scalar (>= 1).');
    end
    
    % Check n_a
    if ~isnumeric(n_a) || ~isscalar(n_a) || floor(n_a) ~= n_a || n_a < 0
        error('n_a must be a non-negative integer scalar (>= 0).');
    end
    
    % Check n_b
    if ~isnumeric(n_b) || ~isscalar(n_b) || floor(n_b) ~= n_b || n_b < 0
        error('n_b must be a non-negative integer scalar (>= 0).');
    end

    

    n_a_states = n*n_a;
    n_b_states = n*n_b;
    n_u_d_states = n;

    n_states = n_a_states + n_b_states + n_u_d_states;

    b_state_indices = n_a_states+(1:n_b_states);

    min_max_range = nan(n_states,2);
    min_max_range(b_state_indices,:) = repmat([0 1],n*n_b,1);
end