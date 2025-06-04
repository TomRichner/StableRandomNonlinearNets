function [r, p] = compute_dependent_variables(a, b, u_d, n_a, n_b, c_SFA)
    % compute r and p from a, b, u_d
    if n_a > 0
        r = relu(u_d - c_SFA .* squeeze_dim(sum(a,2),2));        % Hz, spike rate
    else
        r = relu(u_d);                            % no SFA
    end

    if n_b > 0
        p = r .* squeeze_dim(prod(b,2),2);                       % axonal output
    else
        p = r;                                    % no STD
    end
end