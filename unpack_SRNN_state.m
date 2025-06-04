function [a, b, u_d] = unpack_SRNN_states(X, n, n_a, n_b)
    % Unpack SRNN state vector into named variables
    idx_a_end = n * n_a;

    % ---------- SFA (a–states) ----------
    if n_a > 0
        a = reshape(X(1:idx_a_end), n, n_a);      % n × n_a
    else
        a = [];                                   % no SFA states
    end

    % ---------- STD (b–states) ----------
    idx_b_start = idx_a_end + 1;
    idx_b_end   = idx_a_end + n * n_b;

    if n_b > 0
        b = reshape(X(idx_b_start:idx_b_end), n, n_b); % n × n_b
    else
        b = [];                                   % no STD states
    end

    % ---------- Dendrite (u_d) ----------
    idx_ud_start = idx_b_end + 1;
    idx_ud_end   = idx_b_end + n;
    u_d = X(idx_ud_start:idx_ud_end);
end