function [min_max_range] = get_minMaxRange(params)
    % get_minMaxRange - Defines min/max bounds for state variables based on params struct.
    %
    % Inputs:
    %   params - Struct containing all relevant network and state parameters:
    %            n_E, n_I, n_a_E, n_a_I, n_b_E, n_b_I, n
    %
    % Outputs:
    %   min_max_range - N_sys_eqs x 2 matrix, where N_sys_eqs is the total
    %                   number of state variables. Column 1 is min, Col 2 is max.
    %                   NaN indicates no explicit bound.

    % Extract necessary parameters
    n_E   = params.n_E;
    n_I   = params.n_I;
    n_a_E = params.n_a_E;
    n_a_I = params.n_a_I;
    n_b_E = params.n_b_E;
    n_b_I = params.n_b_I;
    n     = params.n; % Total number of neurons for u_d states

    % Validate inputs (basic checks, can be expanded)
    if ~isstruct(params)
        error('Input must be a params struct.');
    end
    fields_to_check = {'n_E', 'n_I', 'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I', 'n'};
    for i = 1:length(fields_to_check)
        if ~isfield(params, fields_to_check{i}) || ~isnumeric(params.(fields_to_check{i})) || ~isscalar(params.(fields_to_check{i})) || params.(fields_to_check{i}) < 0
            error('Parameter %s is missing, not a non-negative scalar, or invalid.', fields_to_check{i});
        end
        if floor(params.(fields_to_check{i})) ~= params.(fields_to_check{i})
             error('Parameter %s must be an integer.', fields_to_check{i});
        end
    end
    if params.n < (params.n_E + params.n_I)
        error('Total neurons n cannot be less than n_E + n_I.');
    end


    % Calculate the number of states for each block
    num_a_E_states = 0;
    if n_E > 0 && n_a_E > 0
        num_a_E_states = n_E * n_a_E;
    end

    num_a_I_states = 0;
    if n_I > 0 && n_a_I > 0
        num_a_I_states = n_I * n_a_I;
    end

    num_b_E_states = 0;
    if n_E > 0 && n_b_E > 0
        num_b_E_states = n_E * n_b_E;
    end

    num_b_I_states = 0;
    if n_I > 0 && n_b_I > 0
        num_b_I_states = n_I * n_b_I;
    end
    
    num_u_d_states = n; % u_d states are for all n neurons

    % Total number of system equations
    N_sys_eqs = num_a_E_states + num_a_I_states + ...
                num_b_E_states + num_b_I_states + ...
                num_u_d_states;

    if N_sys_eqs == 0 && n > 0 % Case where there might be neurons but no dynamic states (e.g. only u_d if n_a/n_b are all zero)
         min_max_range = nan(n,2); % u_d states default to NaN bounds
         return;
    elseif N_sys_eqs == 0 && n == 0
        min_max_range = []; % No states, no bounds
        return;
    end


    min_max_range = nan(N_sys_eqs, 2); % Initialize with NaN (no bounds)

    current_idx = 0;

    % Bounds for a_E states (typically unbounded: NaN)
    % current_idx remains current_idx + num_a_E_states;

    % Bounds for a_I states (typically unbounded: NaN)
    current_idx = current_idx + num_a_E_states;
    % current_idx remains current_idx + num_a_I_states;
    
    % Bounds for b_E states ([0 1])
    start_b_E = current_idx + num_a_I_states + 1;
    end_b_E   = start_b_E + num_b_E_states - 1;
    if num_b_E_states > 0
        min_max_range(start_b_E:end_b_E, :) = repmat([0 1], num_b_E_states, 1);
    end
    
    % Bounds for b_I states ([0 1])
    start_b_I = end_b_E + 1;
    end_b_I   = start_b_I + num_b_I_states - 1;
    if num_b_I_states > 0
        min_max_range(start_b_I:end_b_I, :) = repmat([0 1], num_b_I_states, 1);
    end

    % Bounds for u_d states (typically unbounded: NaN)
    % The remaining states are u_d, already initialized to NaN.
end