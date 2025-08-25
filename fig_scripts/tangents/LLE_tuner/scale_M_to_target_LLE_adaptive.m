function [scale, LLE_hit, r0, M_adj, DC_adj, c_SFA_adj, adjust_log] = ...
    scale_M_to_target_LLE_adaptive(M, target_LLE, DC, n, n_E, n_I, E_indices, I_indices, ...
                                   n_a_E, tau_a_E, c_SFA, n_b_E, tau_b_E, F_STD, tau_STD, tau_d, varargin)
% Adaptive scaling to hit a target LLE by adjusting DC, E/I block gains, and c_SFA if needed.
%
% Inputs (superset of the original function):
%   M               – connectivity (n×n)
%   target_LLE      – desired Λ_max [1/s]
%   DC              – constant external input level
%   n, n_E, n_I     – sizes
%   E_indices       – indices of excitatory neurons
%   I_indices       – indices of inhibitory neurons
%   n_a_E, tau_a_E, c_SFA, n_b_E, tau_b_E, F_STD, tau_STD, tau_d – model params
%
% Optional name/value:
%   'Tol'           – LLE tolerance (default 1e-4)
%   'MaxIter'       – bisection iterations (default 30)
%   'Verbose'       – logical (default true)
%   'MaxRounds'     – max adjustment rounds (default 5)
%   'MaxExpand'     – max doubling while bracketing scale (default 30)
%   'DCFactor'      – multiplicative DC step (default 1.5)
%   'EEFactor'      – step for E->E increase (default 1.2)
%   'IEMagFactor'   – step to reduce |I->E| (default 0.9)    (moves negative block toward 0)
%   'EIFactor'      – step to reduce E->I (default 0.9)
%   'IIFactor'      – step to reduce |I->I| (default 0.9)
%   'SFAFactor'     – multiplicative c_SFA step (default 0.8)
%
% Outputs:
%   scale           – scalar s such that M_scaled = s*M_adj reaches target Λ_max
%   LLE_hit         – achieved Λ_max at the returned scale
%   r0              – fixed point rates at the returned (M_adj, DC_adj)
%   M_adj           – block-adjusted M used for scaling
%   DC_adj          – adjusted DC
%   c_SFA_adj       – adjusted c_SFA (n×1)
%   adjust_log      – struct with per-round diagnostics/actions

    p = inputParser;
    p.addParameter('Tol', 1e-4);
    p.addParameter('MaxIter', 30);
    p.addParameter('Verbose', true);
    p.addParameter('MaxRounds', 5);
    p.addParameter('MaxExpand', 30);
    p.addParameter('DCFactor', 1.5);
    p.addParameter('EEFactor', 1.2);
    p.addParameter('IEMagFactor', 0.9);
    p.addParameter('EIFactor', 0.9);
    p.addParameter('IIFactor', 0.9);
    p.addParameter('SFAFactor', 0.8);
    p.parse(varargin{:});
    Tol        = p.Results.Tol;
    MaxIter    = p.Results.MaxIter;
    Verbose    = p.Results.Verbose;
    MaxRounds  = p.Results.MaxRounds;
    MaxExpand  = p.Results.MaxExpand;
    DCFactor   = p.Results.DCFactor;
    EEFactor   = p.Results.EEFactor;
    IEMagFact  = p.Results.IEMagFactor;
    EIFact     = p.Results.EIFactor;
    IIFact     = p.Results.IIFactor;
    SFAFact    = p.Results.SFAFactor;

    % Working copies
    M_adj     = M;
    DC_adj    = DC;
    c_SFA_adj = c_SFA(:);

    adjust_log.rounds = struct('round', {}, 'mu_active', {}, 'nE_active', {}, 'nI_active', {}, 'max_u_d0', {}, 'max_r0', {}, 'action', {});
    adjust_log.success = false;

    for round_idx = 0:MaxRounds
        % 1) Try bracketing + bisection on current (M_adj, DC_adj, c_SFA_adj)
        [success, scale_try, LLE_try, r0_try] = try_scale(M_adj, DC_adj, c_SFA_adj);
        if success
            scale   = scale_try;
            LLE_hit = LLE_try;
            r0      = r0_try;
            adjust_log.success = true;
            adjust_log.rounds(end+1).action = 'scale_only';
            if Verbose
                fprintf('[Adaptive] Success with scaling alone at round %d: scale=%.4g, LLE=%.4g\n', round_idx, scale, LLE_hit);
            end
            return
        end

        % 2) Diagnose active set and masked spectral right-edge
        [r0_diag, ~, ~, ~, u_d0_diag] = LLE_analytic_SRNN_robust_extra_stable_fcn( ...
            n, n_E, n_I, M_adj, DC_adj, n_a_E, tau_a_E, c_SFA_adj, ...
            n_b_E, tau_b_E, F_STD, tau_STD, tau_d);
        H = double(r0_diag > 0); % active set
        mu_active = max(real(eig(M_adj * diag(H)))); %#ok<MINV>
        nE_active = sum(H(E_indices) > 0);
        nI_active = sum(H(I_indices) > 0);

        log_entry = struct();
        log_entry.round      = round_idx;
        log_entry.mu_active  = mu_active;
        log_entry.nE_active  = nE_active;
        log_entry.nI_active  = nI_active;
        log_entry.max_u_d0   = max(u_d0_diag);
        log_entry.max_r0     = max(r0_diag);
        log_entry.action     = '';

        if Verbose
            fprintf('[Adaptive] round %d: mu_active=%.4g, nE_active=%d, nI_active=%d, max r0=%.4g\n', ...
                round_idx, mu_active, nE_active, nI_active, max(r0_diag));
        end

        if round_idx == MaxRounds
            % Give up; return best we have (no positive bracket)
            scale   = 1.0;
            LLE_hit = LLE_try;
            r0      = r0_diag;
            adjust_log.rounds(end+1) = log_entry;
            return
        end

        % 3) Adjustment strategy
        if nE_active == 0
            % No excitatory activity → raise DC to recruit E
            DC_adj = DCFactor * DC_adj;
            log_entry.action = sprintf('increase_DC_x%.3g', DCFactor);
        elseif mu_active <= 0
            % Active set lacks a positive-gain mode → tilt E/I
            M_adj(E_indices, E_indices) = EEFactor * M_adj(E_indices, E_indices);   % E->E
            M_adj(E_indices, I_indices) = IEMagFact * M_adj(E_indices, I_indices);  % reduce |I->E|
            M_adj(I_indices, E_indices) = EIFact   * M_adj(I_indices, E_indices);   % reduce E->I
            M_adj(I_indices, I_indices) = IIFact   * M_adj(I_indices, I_indices);   % reduce |I->I|
            log_entry.action = sprintf('EE_x%.2g, IE_mag_x%.2g, EI_x%.2g, II_mag_x%.2g', ...
                EEFactor, IEMagFact, EIFact, IIFact);
        else
            % There is a positive μ over the active set, but still stable → relax SFA
            c_SFA_adj = SFAFact * c_SFA_adj;
            log_entry.action = sprintf('reduce_c_SFA_x%.3g', SFAFact);
        end

        adjust_log.rounds(end+1) = log_entry;
    end

    % Fallback (should not reach)
    scale   = 1.0;
    LLE_hit = NaN;
    r0      = zeros(n,1);

    % ----------------- nested helpers -----------------
    function [ok, scale_out, LLE_out, r0_out] = try_scale(Mloc, DCloc, cSFA_loc)
        lle = @(s) LLE_of_scale(s, Mloc, DCloc, cSFA_loc);

        s_low = 0; s_high = 1;
        L_low = lle(s_low);
        L_high = lle(s_high);

        % Expand upwards
        guard = 0;
        while L_high < target_LLE && guard < MaxExpand
            s_high = 2*s_high + (s_high==0);
            L_high = lle(s_high);
            guard = guard + 1;
        end

        % If baseline already above target, shrink to bracket
        guard = 0;
        while L_low > target_LLE && guard < MaxExpand
            s_high = s_low;
            L_high = L_low;
            s_low  = max(0, 0.5*s_low);
            L_low  = lle(s_low);
            guard  = guard + 1;
        end

        if ~(L_low <= target_LLE && L_high >= target_LLE)
            ok = false; scale_out = NaN; LLE_out = L_high; r0_out = zeros(n,1);
            return
        end

        % Bisection
        for k = 1:MaxIter
            s_mid = 0.5*(s_low + s_high);
            L_mid = lle(s_mid);
            if abs(L_mid - target_LLE) < Tol
                scale_out = s_mid; LLE_out = L_mid;
                [r0_out, ~] = LLE_analytic_SRNN_robust_extra_stable_fcn(n, n_E, n_I, scale_out*Mloc, DCloc, ...
                    n_a_E, tau_a_E, cSFA_loc, n_b_E, tau_b_E, F_STD, tau_STD, tau_d);
                ok = true; return
            end
            if L_mid < target_LLE
                s_low = s_mid; L_low = L_mid;
            else
                s_high = s_mid; L_high = L_mid;
            end
        end

        scale_out = 0.5*(s_low + s_high);
        LLE_out   = lle(scale_out);
        [r0_out, ~] = LLE_analytic_SRNN_robust_extra_stable_fcn(n, n_E, n_I, scale_out*Mloc, DCloc, ...
            n_a_E, tau_a_E, cSFA_loc, n_b_E, tau_b_E, F_STD, tau_STD, tau_d);
        ok = true;
    end

    function L = LLE_of_scale(s, Mloc, DCloc, cSFA_loc)
        [~, L] = LLE_analytic_SRNN_robust_extra_stable_fcn(n, n_E, n_I, s*Mloc, DCloc, ...
            n_a_E, tau_a_E, cSFA_loc, n_b_E, tau_b_E, F_STD, tau_STD, tau_d);
    end
end