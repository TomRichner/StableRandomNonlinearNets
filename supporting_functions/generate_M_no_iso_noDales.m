function [A, EI_vec] = generate_M_no_iso_noDales(n,w,sparsity, EI)
% generates a sparse connection matrix that is strongly connected (all nodes have a directed path to all nodes from graph theory)

    % --------------------------- E / I labels ------------------------------
    EI_vec            = -ones(1,n);          % −1 : inhibitory
    EI_vec(1:round(EI*n)) =  1;              % +1 : excitatory

    % ----------------------- seed-then-thin mask ---------------------------
    E0     = n*(n-1);                        % off-diagonal capacity
    E_keep = round((1-sparsity) * E0);       % total edges to keep

    % feasibility check ----------------------------------------------------
    if E_keep < n
        error('generate_M_no_iso:ImpossibleSparsity', ...
              ['Requested sparsity too high – every neuron must keep at least ' ...
               'one in- and one out-connection, and the network must stay ' ...
               'strongly connected.']);
    end

    mask = false(n,n);                       % true = edge kept
                                             % false = edge removed

    % 0) strongly-connected skeleton : a random Hamiltonian cycle ----------
    perm  = randperm(n);                     % random ordering of neurons
    next  = [perm(2:end) perm(1)];           % successor in the cycle
    mask( sub2ind([n n], perm, next) ) = true; % perm(i) -> perm(i+1)

    % 1) randomly add edges until exact density is reached -----------------
    E_add = E_keep - nnz(mask);              % how many more we still need
    if E_add > 0
        available = find(~mask & ~eye(n));   % remaining off-diag slots
        pick      = available( randperm(numel(available), E_add) );
        mask(pick) = true;
    end
    % ----------------------------------------------------------------------

    % --------------------------- draw weights -----------------------------
    A = randn(n,n);               % Gaussian weights (σ = 1)
    A(~mask) = 0;                 % impose sparsity mask
    A(eye(n,'logical')) = 0;      % clear self-connections for now

    % Dale's law
    A(:, EI_vec== 1) =  abs(A(:, EI_vec== 1));   % outputs from excitatory cells
    A(:, EI_vec==-1) = -abs(A(:, EI_vec==-1));   % outputs from inhibitory cells

    % class-specific scaling
    A(EI_vec==-1, EI_vec== 1) = w.EI .* A(EI_vec==-1, EI_vec== 1);
    A(EI_vec== 1, EI_vec==-1) = w.IE .* A(EI_vec== 1, EI_vec==-1);
    A(EI_vec== 1, EI_vec== 1) = w.EE .* A(EI_vec== 1, EI_vec== 1);
    A(EI_vec==-1, EI_vec==-1) = w.II .* A(EI_vec==-1, EI_vec==-1);

    % optional self-connections (do not affect the in/out guarantees)
    diag_idx                     = 1:n+1:n^2;
    A(diag_idx)                  = w.selfI;        % default inhibitory value
    A(diag_idx(EI_vec==1))       = w.selfE;        % overwrite excitatory

    % ------------------------ Violate Dale's Law ---------------------------
    % This section permutes the columns of A to violate Dale's Law,
    % ensuring that excitatory neurons can form inhibitory connections and
    % vice-versa, while preserving the overall weight distribution.
    % The permutation is done carefully to not disturb self-connections
    % (diagonal) and the connections ensuring the network is strongly
    % connected (the initial Hamiltonian cycle).

    % Get linear indices of connections to preserve
    H_cycle_indices = sub2ind([n n], perm, next);
    
    % Combine with diagonal indices
    protected_indices = union(diag_idx, H_cycle_indices);
    
    % Create a mask for all non-zero elements
    [row_A, col_A] = find(A);
    all_indices = sub2ind(size(A), row_A, col_A);
    
    % Determine which non-zero elements can be permuted
    permutable_indices = setdiff(all_indices, protected_indices);
    
    % Get the values at these permutable locations
    permutable_values = A(permutable_indices);
    
    % Permute the values and place them back into A
    A(permutable_indices) = permutable_values(randperm(numel(permutable_values)));
    % -----------------------------------------------------------------------

    % Sanity check (optional debugging aid)
    G = digraph(mask);          % MATLAB's directed graph object
    comps = conncomp(G, 'Type', 'strong');
    assert(numel(unique(comps)) == 1, 'Graph is not strongly connected');
end

