function J = SRNN_Jacobian(t, X, params)
% SRNN_Jacobian - Analytic Jacobian for the SRNN ODE system with E/I split.
%
% State ordering (consistent with unpack_SRNN_state and SRNN.m):
%   X = [ a_E(:) ; a_I(:) ; b_E(:) ; b_I(:) ; u_d(:) ]
%
% Inputs
%   t       – time (unused, but kept for signature compatibility)
%   X       – current state vector
%   params  – structure produced by package_params.m
%
% Output
%   J       – Jacobian matrix  (N_sys_eqs × N_sys_eqs)

% ---------- Parameters --------------------------------------------------
    n         = params.n;
    M         = params.M;
    tau_d     = params.tau_d;
    c_SFA_all = params.c_SFA; % n x 1
    F_STD_all = params.F_STD; % n x 1
    tau_STD   = params.tau_STD;

    E_indices = params.E_indices;
    I_indices = params.I_indices;
    n_E = params.n_E;
    n_I = params.n_I;

    n_a_E = params.n_a_E; tau_a_E = params.tau_a_E; % tau_a_E is 1 x n_a_E
    n_a_I = params.n_a_I; tau_a_I = params.tau_a_I; % tau_a_I is 1 x n_a_I
    n_b_E = params.n_b_E; tau_b_E = params.tau_b_E; % tau_b_E is 1 x n_b_E
    n_b_I = params.n_b_I; tau_b_I = params.tau_b_I; % tau_b_I is 1 x n_b_I
    
    eps_b_div = 1e-12; % Small number for safe division

% ---------- Unpack State ------------------------------------------------
    [a_E, a_I, b_E, b_I, u_d] = unpack_SRNN_state(X, params);
    % a_E: n_E x n_a_E (or empty)
    % a_I: n_I x n_a_I (or empty)
    % b_E: n_E x n_b_E (or empty)
    % b_I: n_I x n_b_I (or empty)
    % u_d: n x 1

% ---------- Derived Variables (similar to SRNN.m) -----------------------
    u_eff = u_d; % n x 1
    if n_E > 0 && n_a_E > 0 && ~isempty(a_E)
        u_eff(E_indices) = u_eff(E_indices) - c_SFA_all(E_indices) .* sum(a_E, 2);
    end
    if n_I > 0 && n_a_I > 0 && ~isempty(a_I)
        u_eff(I_indices) = u_eff(I_indices) - c_SFA_all(I_indices) .* sum(a_I, 2);
    end
    
    H_all = double(u_eff > 0); % n x 1, derivative of ReLU
    r_all = H_all .* u_eff;    % n x 1, effectively r_all = max(0, u_eff)
                               % but H_all is what's needed for Jacobian

    p_all = r_all; % n x 1
    B_E_prod_per_neuron = ones(n_E, 1); % Product of b_E for each E neuron
    if n_E > 0 && n_b_E > 0 && ~isempty(b_E)
        B_E_prod_per_neuron = prod(b_E, 2); % n_E x 1
        p_all(E_indices) = p_all(E_indices) .* B_E_prod_per_neuron;
    end
    
    B_I_prod_per_neuron = ones(n_I, 1); % Product of b_I for each I neuron
    if n_I > 0 && n_b_I > 0 && ~isempty(b_I)
        B_I_prod_per_neuron = prod(b_I, 2); % n_I x 1
        p_all(I_indices) = p_all(I_indices) .* B_I_prod_per_neuron;
    end

% ---------- Indexing Helpers -------------------------------------------
    num_a_E_states = n_E * n_a_E;
    num_a_I_states = n_I * n_a_I;
    num_b_E_states = n_E * n_b_E;
    num_b_I_states = n_I * n_b_I;
    num_u_d_states = n;

    N_sys_eqs = num_a_E_states + num_a_I_states + num_b_E_states + num_b_I_states + num_u_d_states;

    % Start indices for each block of state variables
    start_a_E = 1;
    start_a_I = start_a_E + num_a_E_states;
    start_b_E = start_a_I + num_a_I_states;
    start_b_I = start_b_E + num_b_E_states;
    start_u_d = start_b_I + num_b_I_states;

% ---------- Allocate Jacobian ------------------------------------------
    J = zeros(N_sys_eqs, N_sys_eqs);

% -----------------------------------------------------------------------
% Block 1: Derivatives of da_E/dt & da_I/dt 
% da_E_dt(e_row, k_col_aE) = (r_all(E_indices(e_row)) - a_E(e_row, k_col_aE)) / tau_a_E(k_col_aE)
% da_I_dt(i_row, k_col_aI) = (r_all(I_indices(i_row)) - a_I(i_row, k_col_aI)) / tau_a_I(k_col_aI)
% -----------------------------------------------------------------------
    % --- 1.1 Derivatives for da_E/dt ---
    if n_E > 0 && n_a_E > 0
        for e_neuron_idx = 1:n_E % Loop over E neurons (local E index)
            neuron_global_idx_E = E_indices(e_neuron_idx);
            if c_SFA_all(neuron_global_idx_E) == 0 % If SFA is off for this E neuron, its a_E derivs are 0
                for k_aE_idx = 1:n_a_E
                    row_J = start_a_E + (e_neuron_idx-1)*n_a_E + k_aE_idx - 1;
                    J(row_J, :) = 0;
                end
                continue; 
            end

            for k_aE_idx = 1:n_a_E % Loop over SFA timescales for this E neuron
                row_J = start_a_E + (e_neuron_idx-1)*n_a_E + k_aE_idx - 1;
                inv_tau_aEk = 1 / tau_a_E(k_aE_idx);
                H_this_E_neuron = H_all(neuron_global_idx_E);
                c_SFA_this_E_neuron = c_SFA_all(neuron_global_idx_E);

                % d(da_E/dt) / da_E'(e',k') : only if e == e'
                for k_aE_prime_idx = 1:n_a_E
                    col_J_aE = start_a_E + (e_neuron_idx-1)*n_a_E + k_aE_prime_idx - 1;
                    delta_k = double(k_aE_idx == k_aE_prime_idx);
                    J(row_J, col_J_aE) = (-c_SFA_this_E_neuron * H_this_E_neuron - delta_k) * inv_tau_aEk;
                end
                % d(da_E/dt) / da_I : is 0 because r_E is not directly func of a_I

                % d(da_E/dt) / du_d(j_global) : only if j_global == neuron_global_idx_E
                col_J_ud_for_this_E = start_u_d + neuron_global_idx_E - 1;
                J(row_J, col_J_ud_for_this_E) = H_this_E_neuron * inv_tau_aEk;
            end
        end
    end

    % --- 1.2 Derivatives for da_I/dt ---
    if n_I > 0 && n_a_I > 0
        for i_neuron_idx = 1:n_I % Loop over I neurons (local I index)
            neuron_global_idx_I = I_indices(i_neuron_idx);
            if c_SFA_all(neuron_global_idx_I) == 0 % If SFA is off for this I neuron
                for k_aI_idx = 1:n_a_I
                    row_J = start_a_I + (i_neuron_idx-1)*n_a_I + k_aI_idx - 1;
                    J(row_J, :) = 0;
                end
                continue;
            end

            for k_aI_idx = 1:n_a_I % Loop over SFA timescales for this I neuron
                row_J = start_a_I + (i_neuron_idx-1)*n_a_I + k_aI_idx - 1;
                inv_tau_aIk = 1 / tau_a_I(k_aI_idx);
                H_this_I_neuron = H_all(neuron_global_idx_I);
                c_SFA_this_I_neuron = c_SFA_all(neuron_global_idx_I);

                % d(da_I/dt) / da_I'(i',k') : only if i == i'
                for k_aI_prime_idx = 1:n_a_I
                    col_J_aI = start_a_I + (i_neuron_idx-1)*n_a_I + k_aI_prime_idx - 1;
                    delta_k = double(k_aI_idx == k_aI_prime_idx);
                    J(row_J, col_J_aI) = (-c_SFA_this_I_neuron * H_this_I_neuron - delta_k) * inv_tau_aIk;
                end
                % d(da_I/dt) / da_E : is 0

                % d(da_I/dt) / du_d(j_global) : only if j_global == neuron_global_idx_I
                col_J_ud_for_this_I = start_u_d + neuron_global_idx_I - 1;
                J(row_J, col_J_ud_for_this_I) = H_this_I_neuron * inv_tau_aIk;
            end
        end
    end

% -----------------------------------------------------------------------
% Block 2: Derivatives of db_E/dt & db_I/dt
% db_E_dt(e,m) = (1-b_E(e,m))/tau_bE(m) - (F_STD_all(E_idx(e))/tau_STD) * p_all(E_idx(e))
% -----------------------------------------------------------------------
    % --- 2.1 Derivatives for db_E/dt ---
    if n_E > 0 && n_b_E > 0
        for e_neuron_idx = 1:n_E % Loop over E neurons
            neuron_global_idx_E = E_indices(e_neuron_idx);
            if F_STD_all(neuron_global_idx_E) == 0 % If STD is off for this E neuron
                for m_bE_idx = 1:n_b_E
                    row_J = start_b_E + (e_neuron_idx-1)*n_b_E + m_bE_idx - 1;
                    J(row_J, :) = 0;
                end
                continue;
            end
            
            K_E_coeff = -F_STD_all(neuron_global_idx_E) / tau_STD;
            H_this_E_neuron = H_all(neuron_global_idx_E);
            c_SFA_this_E_neuron = c_SFA_all(neuron_global_idx_E);
            r_this_E_neuron = r_all(neuron_global_idx_E);
            B_E_prod_this_neuron = B_E_prod_per_neuron(e_neuron_idx);

            for m_bE_idx = 1:n_b_E % Loop over STD timescales for this E neuron
                row_J = start_b_E + (e_neuron_idx-1)*n_b_E + m_bE_idx - 1;
                inv_tau_bEm = 1 / tau_b_E(m_bE_idx);

                % d(db_E/dt) / da_E(e',k') : only if e == e'
                for k_aE_prime_idx = 1:n_a_E
                    col_J_aE = start_a_E + (e_neuron_idx-1)*n_a_E + k_aE_prime_idx - 1;
                    % d(p_E)/da_E = d(r_E * B_E)/da_E = (dr_E/da_E) * B_E
                    dr_daE = -c_SFA_this_E_neuron * H_this_E_neuron;
                    J(row_J, col_J_aE) = J(row_J, col_J_aE) + K_E_coeff * dr_daE * B_E_prod_this_neuron;
                end
                
                % d(db_E/dt) / db_E(e',m') : only if e == e'
                J(row_J, start_b_E + (e_neuron_idx-1)*n_b_E + m_bE_idx - 1) = J(row_J, start_b_E + (e_neuron_idx-1)*n_b_E + m_bE_idx - 1) - inv_tau_bEm; % Explicit leak term
                for m_bE_prime_idx = 1:n_b_E
                    col_J_bE = start_b_E + (e_neuron_idx-1)*n_b_E + m_bE_prime_idx - 1;
                    % d(p_E)/db_E = d(r_E * B_E)/db_E = r_E * (dB_E/db_E_m')
                    dB_dbE_m_prime = calculate_prod_b_excluding_k(b_E(e_neuron_idx,:), m_bE_prime_idx, B_E_prod_this_neuron, eps_b_div);
                    J(row_J, col_J_bE) = J(row_J, col_J_bE) + K_E_coeff * r_this_E_neuron * dB_dbE_m_prime;
                end

                % d(db_E/dt) / du_d(j_global) : only if j_global == neuron_global_idx_E
                col_J_ud_for_this_E = start_u_d + neuron_global_idx_E - 1;
                % d(p_E)/du_d = d(r_E * B_E)/du_d = (dr_E/du_d) * B_E
                dr_dud = H_this_E_neuron;
                J(row_J, col_J_ud_for_this_E) = J(row_J, col_J_ud_for_this_E) + K_E_coeff * dr_dud * B_E_prod_this_neuron;
            end
        end
    end
    
    % --- 2.2 Derivatives for db_I/dt ---
    if n_I > 0 && n_b_I > 0
        for i_neuron_idx = 1:n_I % Loop over I neurons
            neuron_global_idx_I = I_indices(i_neuron_idx);
             if F_STD_all(neuron_global_idx_I) == 0 % If STD is off for this I neuron
                for m_bI_idx = 1:n_b_I
                    row_J = start_b_I + (i_neuron_idx-1)*n_b_I + m_bI_idx - 1;
                    J(row_J, :) = 0;
                end
                continue;
            end

            K_I_coeff = -F_STD_all(neuron_global_idx_I) / tau_STD;
            H_this_I_neuron = H_all(neuron_global_idx_I);
            c_SFA_this_I_neuron = c_SFA_all(neuron_global_idx_I);
            r_this_I_neuron = r_all(neuron_global_idx_I);
            B_I_prod_this_neuron = B_I_prod_per_neuron(i_neuron_idx);

            for m_bI_idx = 1:n_b_I % Loop over STD timescales for this I neuron
                row_J = start_b_I + (i_neuron_idx-1)*n_b_I + m_bI_idx - 1;
                inv_tau_bIm = 1 / tau_b_I(m_bI_idx);

                % d(db_I/dt) / da_I(i',k') : only if i == i'
                for k_aI_prime_idx = 1:n_a_I
                    col_J_aI = start_a_I + (i_neuron_idx-1)*n_a_I + k_aI_prime_idx - 1;
                    dr_daI = -c_SFA_this_I_neuron * H_this_I_neuron;
                    J(row_J, col_J_aI) = J(row_J, col_J_aI) + K_I_coeff * dr_daI * B_I_prod_this_neuron;
                end
                
                % d(db_I/dt) / db_I(i',m') : only if i == i'
                J(row_J, start_b_I + (i_neuron_idx-1)*n_b_I + m_bI_idx - 1) = J(row_J, start_b_I + (i_neuron_idx-1)*n_b_I + m_bI_idx - 1) - inv_tau_bIm; % Explicit leak
                for m_bI_prime_idx = 1:n_b_I
                    col_J_bI = start_b_I + (i_neuron_idx-1)*n_b_I + m_bI_prime_idx - 1;
                    dB_dbI_m_prime = calculate_prod_b_excluding_k(b_I(i_neuron_idx,:), m_bI_prime_idx, B_I_prod_this_neuron, eps_b_div);
                    J(row_J, col_J_bI) = J(row_J, col_J_bI) + K_I_coeff * r_this_I_neuron * dB_dbI_m_prime;
                end

                % d(db_I/dt) / du_d(j_global) : only if j_global == neuron_global_idx_I
                col_J_ud_for_this_I = start_u_d + neuron_global_idx_I - 1;
                dr_dud = H_this_I_neuron;
                J(row_J, col_J_ud_for_this_I) = J(row_J, col_J_ud_for_this_I) + K_I_coeff * dr_dud * B_I_prod_this_neuron;
            end
        end
    end

% -----------------------------------------------------------------------
% Block 3: Derivatives of du_d/dt
% du_d_dt(j_glob) = (-u_d(j_glob) + u_ex(j_glob) + sum_l M(j_glob,l_glob)*p_all(l_glob)) / tau_d
% -----------------------------------------------------------------------
    inv_tau_d = 1 / tau_d;
    for j_global_idx = 1:n % Row of Jacobian (for u_d(j_global_idx))
        row_J = start_u_d + j_global_idx - 1;
        
        % Explicit leak term d(-u_d(j_global_idx)/tau_d) / du_d(j_global_idx)
        J(row_J, start_u_d + j_global_idx - 1) = J(row_J, start_u_d + j_global_idx - 1) - inv_tau_d;

        for l_global_idx = 1:n % Column of Jacobian (for M(j,l)*p(l) term, where l is the presynaptic neuron)
            
            M_jl_inv_tau_d = M(j_global_idx, l_global_idx) * inv_tau_d;
            if M_jl_inv_tau_d == 0, continue; end % No connection, no contribution

            H_l = H_all(l_global_idx);
            c_SFA_l = c_SFA_all(l_global_idx);
            r_l = r_all(l_global_idx);

            % Find if l_global_idx is E or I, and its local index
            is_E_l = false; e_l_idx = 0;
            is_I_l = false; i_l_idx = 0;
            [is_E_member, e_l_idx_temp] = ismember(l_global_idx, E_indices);
            if is_E_member, is_E_l = true; e_l_idx = e_l_idx_temp; end
            [is_I_member, i_l_idx_temp] = ismember(l_global_idx, I_indices);
            if is_I_member, is_I_l = true; i_l_idx = i_l_idx_temp; end

            B_prod_l = 1; % Effective product of b for neuron l
            if is_E_l && n_E > 0 && n_b_E > 0 && ~isempty(b_E)
                B_prod_l = B_E_prod_per_neuron(e_l_idx);
            elseif is_I_l && n_I > 0 && n_b_I > 0 && ~isempty(b_I)
                B_prod_l = B_I_prod_per_neuron(i_l_idx);
            end

            % d(du_d(j)/dt) / da_E(e',k') where E_indices(e') == l_global_idx
            if is_E_l && n_E > 0 && n_a_E > 0
                for k_aE_prime_idx = 1:n_a_E
                    col_J_aE = start_a_E + (e_l_idx-1)*n_a_E + k_aE_prime_idx - 1;
                    % dp(l)/da_E(l,k') = (dr(l)/da_E(l,k')) * B_prod_l
                    dr_daE = -c_SFA_l * H_l;
                    J(row_J, col_J_aE) = J(row_J, col_J_aE) + M_jl_inv_tau_d * dr_daE * B_prod_l;
                end
            end

            % d(du_d(j)/dt) / da_I(i',k') where I_indices(i') == l_global_idx
            if is_I_l && n_I > 0 && n_a_I > 0
                for k_aI_prime_idx = 1:n_a_I
                    col_J_aI = start_a_I + (i_l_idx-1)*n_a_I + k_aI_prime_idx - 1;
                    dr_daI = -c_SFA_l * H_l;
                    J(row_J, col_J_aI) = J(row_J, col_J_aI) + M_jl_inv_tau_d * dr_daI * B_prod_l;
                end
            end
            
            % d(du_d(j)/dt) / db_E(e',m') where E_indices(e') == l_global_idx
            if is_E_l && n_E > 0 && n_b_E > 0 && ~isempty(b_E)
                for m_bE_prime_idx = 1:n_b_E
                    col_J_bE = start_b_E + (e_l_idx-1)*n_b_E + m_bE_prime_idx - 1;
                    % dp(l)/db_E(l,m') = r(l) * (dB_prod_l / db_E(l,m'))
                    dB_dbE_m_prime = calculate_prod_b_excluding_k(b_E(e_l_idx,:), m_bE_prime_idx, B_prod_l, eps_b_div);
                    J(row_J, col_J_bE) = J(row_J, col_J_bE) + M_jl_inv_tau_d * r_l * dB_dbE_m_prime;
                end
            end

            % d(du_d(j)/dt) / db_I(i',m') where I_indices(i') == l_global_idx
            if is_I_l && n_I > 0 && n_b_I > 0 && ~isempty(b_I)
                for m_bI_prime_idx = 1:n_b_I
                    col_J_bI = start_b_I + (i_l_idx-1)*n_b_I + m_bI_prime_idx - 1;
                    dB_dbI_m_prime = calculate_prod_b_excluding_k(b_I(i_l_idx,:), m_bI_prime_idx, B_prod_l, eps_b_div);
                    J(row_J, col_J_bI) = J(row_J, col_J_bI) + M_jl_inv_tau_d * r_l * dB_dbI_m_prime;
                end
            end
            
            % d(du_d(j)/dt) / du_d(l_global_idx)
            col_J_ud = start_u_d + l_global_idx - 1;
            % dp(l)/du_d(l) = (dr(l)/du_d(l)) * B_prod_l
            dr_dud = H_l;
            J(row_J, col_J_ud) = J(row_J, col_J_ud) + M_jl_inv_tau_d * dr_dud * B_prod_l;
        end
    end
end

function prod_val = calculate_prod_b_excluding_k(b_row_neuron, k_exclude_idx_in_row, B_neuron_val_total, epsilon)
% Calculates prod_{m != k_exclude_idx_in_row} b_row_neuron(m)
% b_row_neuron: 1 x n_b_specific (e.g., 1 x n_b_E) vector of b states for a single neuron
% k_exclude_idx_in_row: index (1 to n_b_specific) of the b_m to exclude
% B_neuron_val_total: precomputed prod(b_row_neuron) for that neuron
% epsilon: small value for checking if b_row_neuron(k_exclude_idx_in_row) is zero

    num_b_timescales_for_neuron = length(b_row_neuron);

    if num_b_timescales_for_neuron == 0 % Should not happen if called
        prod_val = 1; 
        return;
    end
    if num_b_timescales_for_neuron == 1 % If only one b, derivative of product B w.r.t b1 is 1 (empty product)
        prod_val = 1;
        return;
    end

    val_b_k_exclude = b_row_neuron(k_exclude_idx_in_row);
    if abs(val_b_k_exclude) > epsilon
        prod_val = B_neuron_val_total / val_b_k_exclude;
    else
        % val_b_k_exclude is zero or close to it.
        % Calculate product of other elements directly.
        indices = true(1, num_b_timescales_for_neuron);
        indices(k_exclude_idx_in_row) = false;
        prod_val = prod(b_row_neuron(indices));
    end
end 