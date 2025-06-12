function [A, EI_vec] = generate_M(nvars,w,sparseness, EI)

    % EI is now a parameter representing the fraction of neurons that are excitatory
    EI_vec = -1*ones(1,nvars); % vector of E vs I.  E is coded as +1, I is coded as -1, default all to -1
    EI_vec(1:round(EI*nvars)) = 1; % set the first round(EI*nvars) to be excitatory
    A=1*randn(nvars,nvars); % generate a gaussian random connection matrix with standard deviation 1.0 (scaled later)
    A(rand(nvars^2,1)<sparseness) = 0; % delete some connections to make it sparse

    % A(:,1) = A(:,1)*2; % make node 1 drive the network harder, prettier example
    % A(1,:) = A(1,:)/1.5; % make ohter nodes drive node 1 weaker

    A(eye(size(A),'logical')) = 0; % zero main diagonal, no self connections
    A(:,EI_vec==1) = abs(A(:,EI_vec==1)); % enforce excitatory connections from excitatory neurons
    A(:,EI_vec==-1) = -abs(A(:,EI_vec==-1)); % enforce inhibitory connections from inhibitory neurons Dale's law

    A(EI_vec==-1,EI_vec==1) = w.EI.*A(EI_vec==-1,EI_vec==1); % apply connection weights E to I
    A(EI_vec==1,EI_vec==-1) = w.IE.*A(EI_vec==1,EI_vec==-1); % I to E
    A(EI_vec==1,EI_vec==1) = w.EE.*A(EI_vec==1,EI_vec==1); % E to E
    A(EI_vec==-1,EI_vec==-1) = w.II.*A(EI_vec==-1,EI_vec==-1); % I to I

    % max_w = 0.9;
    % A(A>=max_w) = mod(A(A>=max_w),max_w); % no weights over 0.9;

    A(diag(EI_vec==1)) = w.selfE; % could add self connections to excitatory neurons here, turned off

    A(diag(EI_vec==-1)) = w.selfI;

end

