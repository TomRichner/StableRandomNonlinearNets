function J = SRNN_Jacobian(t, X, params)
% Analytic Jacobian for the SRNN ODE system (see SRNN.m)
%
% State ordering (see unpack_SRNN_state.m):
%   X = [ a(:) ; b(:) ; u_d ]
%
% ReLU derivative is treated as a Heaviside step:  H(s) = 1 for s>0, 0 else.
% For the problem sizes of interest ( ≤ 600 states ) a dense matrix is fine.
%
% Inputs
%   t       – time (unused, but kept for signature compatibility)
%   X       – current state vector
%   params  – structure produced by package_params.m
%
% Output
%   J       – Jacobian matrix  (N_sys_eqs × N_sys_eqs)

% ---------- parameters --------------------------------------------------
n        = params.n;
n_a      = params.n_a;
n_b      = params.n_b;

tau_a    = params.tau_a(:)';     % 1 × n_a
tau_b    = params.tau_b(:)';     % 1 × n_b
tau_d    = params.tau_d;

M        = params.M;             % n × n
c_SFA    = params.c_SFA(:);      % n × 1
F_STD    = params.F_STD(:);      % n × 1
tau_STD  = params.tau_STD;

eps_b_div = 1e-12; % Small number to prevent division by zero for b_ik derivatives

% ---------- unpack state ------------------------------------------------
[a, b, u_d] = unpack_SRNN_state(X, n, n_a, n_b);

% ---------- derived variables ------------------------------------------
if n_a > 0
    s_arg = u_d - c_SFA .* sum(a, 2);        % ReLU argument
else
    s_arg = u_d;
end
H = double(s_arg > 0);                       % derivative of ReLU
r = H .* s_arg;                              % spike rate

if n_b > 0
    B = prod(b, 2);                      % ∏_m b_{i,m}
else
    B = ones(n, 1);
end

% ---------- indexing helpers -------------------------------------------
N_sys_eqs = n * (n_a + n_b + 1);

if n_a > 0
    idx_a_all = reshape(1:n*n_a,  n, n_a);           % (i,k) → linear index
else
    idx_a_all = [];
end
if n_b > 0
    idx_b_all = n*n_a + reshape(1:n*n_b, n, n_b);
else
    idx_b_all = [];
end
idx_ud_all = n*n_a + n*n_b + (1:n);                  % dendrite indices

% ---------- allocate Jacobian ------------------------------------------
J = zeros(N_sys_eqs, N_sys_eqs);

% -----------------------------------------------------------------------
% 1) Adaptation states a_{i,k}
% -----------------------------------------------------------------------
if n_a > 0
    for i = 1:n
        if c_SFA(i) == 0
            for k = 1:n_a
                row = idx_a_all(i, k);
                J(row, :) = 0;
            end
            continue;
        end

        for k = 1:n_a
            row = idx_a_all(i, k);
            tauAk = tau_a(k);
            
            % ∂(da_{i,k}/dt)/∂u_{d,i}
            J(row, idx_ud_all(i)) =  H(i) / tauAk;
            
            % ∂(da_{i,k}/dt)/∂a_{i,k′}
            for k2 = 1:n_a
                col = idx_a_all(i, k2);
                delta = double(k == k2);
                J(row, col) = (-delta - c_SFA(i) * H(i)) / tauAk;
            end
        end
    end
end

% -----------------------------------------------------------------------
% 2) Depression states b_{i,m}
% -----------------------------------------------------------------------
if n_b > 0
    for i = 1:n
        if F_STD(i) == 0
            for m = 1:n_b
                row = idx_b_all(i, m);
                J(row, :) = 0;
            end
            continue;
        end

        term2_coeff = -F_STD(i) / tau_STD;

        for m = 1:n_b
            row = idx_b_all(i, m);
            tauBm = tau_b(m);
            
            % explicit leak term
            J(row, idx_b_all(i, m)) = -1 / tauBm;
            
            % derivative via p_i
            % (a) w.r.t. u_{d,i}
            J(row, idx_ud_all(i)) = J(row, idx_ud_all(i)) ...
                + term2_coeff * B(i) * H(i);
            
            % (b) w.r.t. a_{i,k}
            if n_a > 0
                for k = 1:n_a
                    J(row, idx_a_all(i, k)) = J(row, idx_a_all(i, k)) ...
                        + term2_coeff * c_SFA(i) * H(i) * B(i);
                end
            end
            
            % (c) w.r.t. b_{i,k′}
            for k = 1:n_b
                if k == m
                    continue;
                end
                col_b = idx_b_all(i, k);
                denom = max(b(i, k), eps_b_div);
                J(row, col_b) = J(row, col_b) ...
                    + term2_coeff * r(i) * B(i) * calculate_prod_b_excluding_k(b(i,:), k, B(i), eps_b_div) / denom;
            end
        end
    end
end

% -----------------------------------------------------------------------
% 3) Dendritic variables u_{d,i}
% -----------------------------------------------------------------------
for i = 1:n
    row = idx_ud_all(i);
    
    for j = 1:n
        scale = M(i, j) / tau_d;    % common factor
        
        % (a) w.r.t. u_{d,j}
        J(row, idx_ud_all(j)) = J(row, idx_ud_all(j)) + scale * B(j) * H(j);
        
        % (b) w.r.t. a_{j,k}
        if n_a > 0
            for k = 1:n_a
                J(row, idx_a_all(j, k)) = J(row, idx_a_all(j, k)) ...
                    - scale * c_SFA(j) * H(j) * B(j);
            end
        end
        
        % (c) w.r.t. b_{j,m}
        if n_b > 0
            for m = 1:n_b
                denom = max(b(j, m), 1e-12);
                J(row, idx_b_all(j, m)) = J(row, idx_b_all(j, m)) ...
                    + scale * r(j) * B(j) * calculate_prod_b_excluding_k(b(j,:), m, B(j), eps_b_div) / denom;
            end
        end
    end
    
    % explicit leak term:  -u_{d,i}/τ_d
    J(row, idx_ud_all(i)) = J(row, idx_ud_all(i)) - 1 / tau_d;
end
end

function prod_val = calculate_prod_b_excluding_k(b_row_neuron, k_exclude_idx, B_neuron_val, epsilon)
% Calculates prod_{m != k_exclude_idx} b_row_neuron(m)
% b_row_neuron: 1 x n_b vector of b states for a single neuron
% k_exclude_idx: index (1 to n_b) of the b_m to exclude from product
% B_neuron_val: precomputed prod(b_row_neuron)
% epsilon: small value for checking if b_row_neuron(k_exclude_idx) is zero

    n_b_local = length(b_row_neuron);
    if n_b_local == 0
        prod_val = 1; % Or handle as error, but consistent with B=1 if no b states
        return;
    end
    if n_b_local == 1 
        % B = b_1. dB/db_1 = 1.
        % k_exclude_idx must be 1. Product of "others" is 1 (empty product).
        prod_val = 1;
        return;
    end

    if abs(b_row_neuron(k_exclude_idx)) > epsilon
        prod_val = B_neuron_val / b_row_neuron(k_exclude_idx);
    else
        % b_row_neuron(k_exclude_idx) is zero or close to it.
        % Calculate product of other elements directly.
        indices = true(1, n_b_local);
        indices(k_exclude_idx) = false;
        prod_val = prod(b_row_neuron(indices));
    end
end 