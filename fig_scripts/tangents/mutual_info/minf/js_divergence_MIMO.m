function D = js_divergence_MIMO(A, B, L)
    % A, B : n_ch x time
    % L    : number of histogram bins
    % D    : 1×T vector of Jensen–Shannon divergences per frame (bits)
    
    T = size(A,3);
    edges = linspace( double(min(A(:),B(:))), ...
                      double(max(A(:),B(:))), L+1 );
    
    D = zeros(1,T);
    for t = 1:T
        p = histcounts( double(A(:,:,t)), edges, 'Normalization','probability');
        q = histcounts( double(B(:,:,t)), edges, 'Normalization','probability');
    
        m   = 0.5*(p+q);
        % Jensen–Shannon divergence (information‑theoretic, symmetric, finite)
        D(t) = 0.5*sum(  p.*log2((p+eps)./(m+eps)) + ...
                         q.*log2((q+eps)./(m+eps))   );
    end

end