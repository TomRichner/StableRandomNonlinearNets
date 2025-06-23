function [r0_analytic, LLE_analytic] = LLE_analytic_SRNN_fcn(n, nE, nI, M, u_ex_dc, n_a_E, tau_a_E, c_SFA, n_b_E, tau_b_E, F_STD, tau_STD, tau_d)
%% -----------------------------------------------------------
%  File:   LLE_analytic_SRNN_fcn.m
%
%  Analytic fixed‑point and Lyapunov calculations for the
%  SRNN rate model with ReLU, multi‑τ SFA and single‑τ STD.
%
%  Assumptions:
%     • identical SFA / STD parameters for all excitatory neurons
%     • n_b = 1  (single STD time‑constant)
%     • all neurons receive the same constant DC input u_ex
%     • all excitatory neurons are active at the fixed point
%       (r0_E > 0) ; inhibitory neurons may be silent or active
%
%  -----------------------------------------------------------

% For analytic solution, we require single timescale STD on E neurons
if n_b_E ~= 1
    warning('Analytic solution requires n_b_E = 1. Skipping calculation.');
    r0_analytic = nan(n,1);
    LLE_analytic = NaN;
    return;
end

% External drive (DC) -----------------------------------------------------
u_ex = u_ex_dc * ones(n,1);   % identical DC to all neurons

tol_fp = 1e-9;             % fixed‑point tolerance
maxIt  = 10000;

%% ------------------ 2.   FIXED‑POINT ITERATION -------------------------

% helper: gamma = F_STD * tau_b / tau_STD (only non‑zero on E)
gamma = F_STD .* (tau_b_E/tau_STD);

alpha  = 0.1;          % 0<alpha≤1  (smaller ⇒ safer, slower)
r      = max(0,u_ex);  % initial guess
for it = 1:maxIt
    p = r ./ (1 + gamma .* r);
    u_d = u_ex + M*p;
    u_eff = u_d - c_SFA .* n_a_E .* r;
    r_raw = max(0, u_eff);
    r_new = (1-alpha)*r + alpha*r_raw;   % Damped update
    
    if norm(r_new - r,inf) < tol_fp
        r = r_new;  break
    end
    
    if it==maxIt,  error('Fixed point did not converge'); end

    r = r_new;
end
% fprintf('Analytic Fixed‑point solved in %d iterations\n',it);

r0 = r;
r0_analytic = r0;
p0 = r0 ./ (1 + gamma .* r0);
h0 = 1 ./ (1 + gamma .* r0);
b0 = h0;                     % single‑τ => b0 == h0

% sanity: print max firing rate
% fprintf('max(r0_analytic) = %.4g Hz\n',max(r0));

%% ------------------ 3.   BUILD POLYNOMIAL COEFFICIENTS -----------------

% Convenience: same parameters for all *active* E neurons
E_indices = 1:nE;
I_indices = (nE+1):n;

active_E_indices = find(r0(E_indices) > 0);
if isempty(active_E_indices)
    warning('No excitatory neurons are active at the fixed point. Analytic LLE calculation may not be valid.');
    LLE_analytic = NaN;
    return;
end
active_E_indices = E_indices(active_E_indices); % global indices

first_active_E = active_E_indices(1);

if any(c_SFA(active_E_indices)~=c_SFA(first_active_E)) || ...
   any(gamma(active_E_indices)~=gamma(first_active_E))
    warning('Heterogeneous adaptation parameters – analytic block‑diagonal derivation is needed.');
end
cS  = c_SFA(first_active_E);
rS  = r0(first_active_E);
bS  = b0(first_active_E);
hS  = h0(first_active_E);
% K_val = F_STD(first_active_E)*rS/(bS*tau_STD);
% gammaS = gamma(first_active_E);

% ---- A(λ)   (SFA) -------------------------------------------------------
DenA = 1;                       % coefficients, highest power first
for k = 1:n_a_E
    DenA = conv(DenA,[tau_a_E(k) 1]);     % multiply (τ_a λ + 1)
end

NumerA = DenA;                  % start with "1" -> becomes (tau_a*lambda+1)
A_numer_sum_term = 0;
for k = 1:n_a_E
    % divide by (τ_a λ + 1) to get prod_{j≠k}
    Partial = deconv(DenA,[tau_a_E(k) 1]);
    A_numer_sum_term = polyadd(A_numer_sum_term, cS*Partial);   %#ok<*NOPRT>
end
NumerA = polyadd(NumerA, A_numer_sum_term);

% ---- B(λ)   (STD) -------------------------------------------------------
K   = F_STD(first_active_E)*rS/(bS*tau_STD);   % dimensionless
DenB   = [tau_b_E 1];              % τ_b*λ + 1
NumerB = [tau_b_E 1+K*tau_b_E];    % τ_b*λ + 1 + K*τ_b

% ---- C(λ) = (λ τ_d + 1) -----------------------------------------------
Cpoly = [tau_d 1];

% helper for poly multiplication
ABC_num = conv(conv(Cpoly,NumerA),NumerB);
ABD_den = conv(conv(DenA,DenB), Cpoly);

%% ------------------ 4.   EIGEN‑LOOP  &  LLE ----------------------------
eigM = eig(M);           % μ_k
LLE  = -Inf;

% The characteristic equation is P(lambda) = 0 for each mu_k.
% P(lambda) = ABD_den(lambda) - mu_k * hS * ABC_num(lambda)
% Note: The transfer function from p to u_d is M/(tau_d*lambda+1).
% The derivation in the original script seemed to have combined polynomials
% in a way that led to a different characteristic equation. This has been
% corrected to reflect the system dynamics more accurately.
% The logic for constructing the equation is:
% delta_r = G(lambda) * M * delta_p
% delta_p = H(lambda) * delta_r
% --> 1 = G(lambda) * mu_k * H(lambda) where mu_k is eigenvalue of M
% G(lambda) = 1 / ( (tau_d*lambda+1) * A(lambda) )
% H(lambda) = hS * B(lambda)
% --> (tau_d*lambda+1)*A(lambda) = mu_k * hS * B(lambda)
% --> (tau_d*lambda+1)*NumerA/DenA = mu_k * hS * NumerB/DenB
% --> DenB * Cpoly * NumerA = mu_k * hS * NumerB * DenA

DenA_DenB = conv(DenA, DenB);
LHS_poly = conv(Cpoly, DenA_DenB);

for k = 1:n
    mu = eigM(k);
    
    RHS_poly = mu * hS * conv(NumerA, NumerB);

    len = max(numel(LHS_poly),numel(RHS_poly));
    p_left  = [zeros(1,len-numel(LHS_poly))  LHS_poly];
    p_right = [zeros(1,len-numel(RHS_poly))  RHS_poly];

    polyEq  = p_left - p_right;
    
    lambdaRoots = roots(polyEq);
    LLE = max(LLE, max(real(lambdaRoots)));
end
% fprintf('Analytic largest Lyapunov exponent  Λ_max = %+8.5f  1/s\n',LLE);
LLE_analytic = LLE;

end


%% ------------------ 5.   HELPER (polyadd) ------------------------------
function s = polyadd(a,b)
% add two polynomial coefficient vectors (highest power first)
if numel(a)<numel(b)
    a = [zeros(1,numel(b)-numel(a)) a];
elseif numel(b)<numel(a)
    b = [zeros(1,numel(a)-numel(b)) b];
end
s = a+b;
end
