function [scale, LLE_hit, r0] = scale_M_to_target_LLE_noSTD(M, target_LLE, DC, n, n_E, n_I, n_a_E, tau_a_E, c_SFA, n_b_E, tau_b_E, F_STD, tau_STD, tau_d, varargin)
% General scaling to hit a target LLE using the robust analytic solver.
% Supports SFA, STD, both or none (as set by n_a_E, n_b_E, c_SFA, F_STD).
%
% Inputs
%   M           – connectivity (n×n)
%   target_LLE  – desired Λ_max [1/s]
%   DC          – constant external input level
%   n, n_E, n_I – sizes (n_E+n_I must equal n)
%   n_a_E       – # of SFA timescales on E neurons (0 ⇒ no SFA)
%   tau_a_E     – 1×n_a_E SFA time constants (empty if n_a_E==0)
%   c_SFA       – n×1 SFA strength per neuron (0 for neurons without SFA)
%   n_b_E       – # of STD timescales on E neurons (0 ⇒ no STD)
%   tau_b_E     – 1×n_b_E STD time constants (empty if n_b_E==0)
%   F_STD       – n×1 STD gain per neuron (0 where STD off)
%   tau_STD     – STD onset time constant (scalar)
%   tau_d       – dendritic time constant (scalar)
%
% Outputs
%   scale       – scalar s such that M_scaled = s*M reaches target Λ_max
%   LLE_hit     – achieved Λ_max at the returned scale
%   r0          – fixed point rates at the returned scale

    p = inputParser;
    p.addParameter('Tol', 1e-4);
    p.addParameter('MaxIter', 30);
    p.parse(varargin{:});
    Tol = p.Results.Tol; MaxIter = p.Results.MaxIter;

    % Handle pure linear no-adaptation/no-STD case analytically
    no_SFA = (n_a_E==0) || all(c_SFA==0);
    no_STD = (n_b_E==0) || all(F_STD==0);

    if no_SFA && no_STD
        mu_max = max(real(eig(M)));
        % Λ_max(s) = (s*Re(μ_max) - 1)/τ_d
        if mu_max <= 0
            scale = 0;
            r0 = zeros(n,1);
            LLE_hit = -1/tau_d;
            return;
        end
        scale = (tau_d*target_LLE + 1) / mu_max;
        % Ensure non-negative scale
        scale = max(0, scale);
        % Report achieved LLE
        LLE_hit = (scale*mu_max - 1)/tau_d;
        % Rough fixed point for ReLU(+DC): solved by robust analytic for consistency
        [r0, ~] = LLE_analytic_SRNN_robust_extra_stable_fcn(n, n_E, n_I, scale*M, DC, n_a_E, tau_a_E, c_SFA, n_b_E, tau_b_E, F_STD, tau_STD, tau_d);
        return;
    end

    % Otherwise use robust analytic LLE to define Λ_max(s)
    lle = @(s) LLE_of_scale(s);

    % Initial bracket
    s_low = 0; s_high = 1;
    L_low = lle(s_low);
    L_high = lle(s_high);

    % Expand upward if needed
    it_guard = 0;
    while L_high < target_LLE && it_guard < 30
        s_high = 2*s_high + (s_high==0);
        L_high = lle(s_high);
        it_guard = it_guard + 1;
    end

    % If baseline already above target, shrink to bracket
    it_guard = 0;
    while L_low > target_LLE && it_guard < 30
        s_high = s_low;          % move high down
        L_high = L_low;
        s_low = max(0, 0.5*s_low);  % shrink low toward 0
        L_low = lle(s_low);
        it_guard = it_guard + 1;
    end

    % Bisection
    for k = 1:MaxIter
        s_mid = 0.5*(s_low + s_high);
        L_mid = lle(s_mid);
        if abs(L_mid - target_LLE) < Tol
            scale = s_mid; LLE_hit = L_mid;
            [r0, ~] = LLE_analytic_SRNN_robust_extra_stable_fcn(n, n_E, n_I, scale*M, DC, n_a_E, tau_a_E, c_SFA, n_b_E, tau_b_E, F_STD, tau_STD, tau_d);
            return;
        end
        if L_mid < target_LLE
            s_low = s_mid; L_low = L_mid;
        else
            s_high = s_mid; L_high = L_mid;
        end
    end

    scale = 0.5*(s_low + s_high);
    LLE_hit = lle(scale);
    [r0, ~] = LLE_analytic_SRNN_robust_extra_stable_fcn(n, n_E, n_I, scale*M, DC, n_a_E, tau_a_E, c_SFA, n_b_E, tau_b_E, F_STD, tau_STD, tau_d);

    function L = LLE_of_scale(s)
        [~, L] = LLE_analytic_SRNN_robust_extra_stable_fcn(n, n_E, n_I, s*M, DC, n_a_E, tau_a_E, c_SFA, n_b_E, tau_b_E, F_STD, tau_STD, tau_d);
    end
end