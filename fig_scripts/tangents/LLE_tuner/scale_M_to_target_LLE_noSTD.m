function [scale, LLE_hit, r0] = scale_M_to_target_LLE_noSTD(M, target_LLE, DC, n_a, tau_a, c_SFA, tau_d, varargin)
% Find scale s so that LLE_analytic(M*s) ~= target_LLE with STD disabled.
% Inputs: M, target_LLE [1/s], DC, SFA params, tau_d
% Outputs: scale, LLE_hit, r0 at the hit.

    p = inputParser;
    p.addParameter('Tol', 1e-4);
    p.addParameter('MaxIter', 30);
    p.parse(varargin{:});
    Tol = p.Results.Tol; MaxIter = p.Results.MaxIter;

    n = size(M,1);
    lle = @(s) LLE_of_scale(s);

    % Quick closed-form guess if no SFA (all c_SFA==0 or n_a==0)
    if (n_a==0) || all(c_SFA==0)
        mu_max = max(real(eig(M)));
        if mu_max <= 0
            scale = 0; [r0, LLE_hit] = LLE_analytic_SRNN_no_STD(n, 0*M, DC, n_a, tau_a, c_SFA, tau_d); return;
        end
        s0 = (tau_d*target_LLE + 1) / mu_max;
        if s0 > 0, s_low = 0; s_high = 2*s0; else, s_low = 0; s_high = 1; end
    else
        s_low = 0; s_high = 1;
    end

    % Ensure bracket: lle(s_low) <= target <= lle(s_high)
    L_low = lle(s_low);
    L_high = lle(s_high);
    it_guard = 0;
    while L_high < target_LLE && it_guard < 30
        s_high = 2*s_high + (s_high==0);
        L_high = lle(s_high);
        it_guard = it_guard + 1;
    end

    % Bisection
    for k = 1:MaxIter
        s_mid = 0.5*(s_low + s_high);
        L_mid = lle(s_mid);
        if abs(L_mid - target_LLE) < Tol
            scale = s_mid; LLE_hit = L_mid; [r0, ~] = LLE_analytic_SRNN_no_STD(n, s_mid*M, DC, n_a, tau_a, c_SFA, tau_d); return;
        end
        if L_mid < target_LLE
            s_low = s_mid; L_low = L_mid;
        else
            s_high = s_mid; L_high = L_mid;
        end
    end
    scale = 0.5*(s_low + s_high); LLE_hit = lle(scale);
    [r0, ~] = LLE_analytic_SRNN_no_STD(n, scale*M, DC, n_a, tau_a, c_SFA, tau_d);

    function L = LLE_of_scale(s)
        [~, L] = LLE_analytic_SRNN_no_STD(n, s*M, DC, n_a, tau_a, c_SFA, tau_d);
    end
end