function test_M_iso()
% Sanity checks for generate_M_no_iso.m
% ----------------------------------------------------------

rng(42);                                 % deterministic tests
nTrials = 100;

for trial = 1:nTrials
    % ------------------- random problem size -----------------
    n          = randi([4 200]);          % 4 ≤ n ≤ 200
    EI_fraction = rand();                 % excitatory share  ∈ (0,1)
    
    % feasible sparsity:   ≤ 1 - 1/(n-1)
    maxSpar    = 1 - 1/(n-1)-1e-9;     % tiny safety margin
    sparsity   = rand() * maxSpar;        % uniform feasible value
    
    % weight scalings (any positive numbers will do for the test)
    w = struct('EI',1,'IE',1,'EE',1,'II',1,'selfE',0,'selfI',0);
    
    % ------------------- call the generator ------------------
    [A,~] = generate_M_no_iso(n,w,sparsity,EI_fraction);
    
    % remove self-loops from consideration
    A(1:n+1:end) = 0;
    
    % 1) each neuron has ≥1 outgoing & incoming edge ------------
    assert(all(sum(A~=0,2) >= 1), 'Out-degree test failed');
    assert(all(sum(A~=0,1) >= 1), 'In-degree  test failed');
    
    % 2) graph is strongly connected ---------------------------
    G      = digraph(A~=0);              % boolean adjacency
    comps  = conncomp(G,'Type','strong');
    assert(numel(unique(comps))==1, 'Graph is not strongly connected');
    
    % 3) realised sparsity matches request ---------------------
    offDiagTotal   = n*(n-1);
    realisedSpars  = (offDiagTotal - nnz(A)) / offDiagTotal;
    assert(abs(realisedSpars - sparsity) <= 1/offDiagTotal, ...
           'Sparsity mismatch');
end

% 4) impossible sparsity should raise an error ----------------
n  = 10;
w  = struct('EI',1,'IE',1,'EE',1,'II',1,'selfE',0,'selfI',0);
try
    generate_M_no_iso(n,w, 0.999, 0.5);  % infeasible: too many zeros
    error('Impossible sparsity did NOT raise an error');
catch ME
    assert(strcmp(ME.identifier, 'generate_M_no_iso:ImpossibleSparsity'), ...
           'Unexpected error identifier for impossible sparsity');
end

disp('All generate_M_no_iso tests passed');
end
