function [r0_analytic, LLE_analytic] = LLE_analytic_SRNN_robust_fcn(n, nE, nI, M, u_ex_dc, n_a, tau_a, c_SFA, n_b, tau_b, F_STD, tau_STD, tau_d)
%% -----------------------------------------------------------
%  File :  LLE_analytic_SRNN_robust_fcn.m
%  Author: ChatGPT (OpenAI o3) & Gemini 2.5 Pro
%  Date :  2025‑06‑12
%
%  – Finds the steady‑state (fixed point) of the SRNN rate model
%    with ReLU, multi‑τ SFA  and  *multi‑τ STD*.
%  – Analytically computes the largest Lyapunov exponent Λ_max
%    at that fixed point.
%
%  The fixed‑point solver is the robust two‑stage algorithm:
%     Stage‑1  damped Picard + Anderson (depth 2)
%     Stage‑2  Jacobian‑free Newton–Krylov (GMRES)
%  -----------------------------------------------------------


%% ---------------- 1.  PARAMETER SETUP ------------------------

% Ensure tau vectors are columns
tau_a = tau_a(:);
tau_b = tau_b(:);

% External drive (constant DC)
u_ex = u_ex_dc * ones(n,1);


%% ---------------- 2.  FIXED‑POINT SEARCH ------------------------------

[r0,b0_mat,h0,infoFP] = fixed_point_multiSTD( ...
        u_ex,M,c_SFA,tau_a,F_STD,tau_b,tau_STD,n_a);

fprintf('FP (robust) | Picard it=%d  Newton it=%d  ‖g‖∞=%.1e   max(r0)=%.4g\n',...
        infoFP.it_picard,infoFP.it_newton,infoFP.res_final,max(r0));

%% ---------------- 3.  ANALYTIC LLE  (multi‑τ STD) --------------------

% Choose a representative *active excitatory* neuron
idxA = (r0>0);
iRef = find(idxA & (c_SFA>0 | F_STD>0),1,'first');
if isempty(iRef)
    warning('No active, adapting excitatory neuron – increase DC or weights. Cannot compute analytic LLE.');
    r0_analytic = r0;
    LLE_analytic = NaN;
    return;
end

% Shorthand for that neuron
rS  = r0(iRef);
bS  = b0_mat(iRef,:).';                 % n_b × 1
hS  = prod(bS);                         % product of bS components
cS  = c_SFA(iRef);
Kvec = F_STD(iRef) * rS * tau_b ./ (bS * tau_STD);  % K_m  (n_b×1)

% ---------- polynomials :   SFA (A)  and  STD (B) ----------------------
DenA = 1;
for k = 1:n_a
    DenA = conv(DenA,[tau_a(k) 1]);            % ∏ (τ_ak λ+1)
end
NumerA = DenA;
if cS > 0 && n_a > 0
    A_numer_sum_term = 0;
    for k = 1:n_a
        Partial = deconv(DenA,[tau_a(k) 1]);
        A_numer_sum_term = polyadd(A_numer_sum_term, cS*Partial);
    end
    NumerA = polyadd(NumerA, A_numer_sum_term);
end


% multi‑τ  STD
DenB = 1;
if n_b > 0
    for m = 1:n_b
        DenB = conv(DenB,[tau_b(m) 1]);          % ∏ (τ_bm λ+1)
    end
    NumerB = DenB;
    for m = 1:n_b
        Partial = deconv(DenB,[tau_b(m) 1]);
        NumerB  = polyadd(NumerB, Kvec(m)*Partial);
    end
else
    NumerB = 1;
end


% dendrite factor  (λ τ_d + 1)
Cpoly  = [tau_d 1];

% Full transfer function polynomials
% G(λ) = 1 / ( (τ_d*λ+1) * A(λ) ) -> (DenA*DenC) / (NumerA*NumerC)
% H(λ) = hS * B(λ) -> hS * NumerB/DenB
% Characteristic eq: 1 = G(λ) * μ * H(λ)
% 1 = (1/ (Cpoly*A(λ))) * mu * hS * B(λ)
% Cpoly(λ)*A(λ) = mu * hS * B(λ)
% Cpoly(λ)*NumerA(λ)/DenA(λ) = mu * hS * NumerB(λ)/DenB(λ)
% Cpoly(λ)*NumerA(λ)*DenB(λ) = mu * hS * NumerB(λ)*DenA(λ)

LHS_poly = conv(conv(Cpoly, NumerA), DenB);
RHS_poly_base = conv(NumerB, DenA);


% ---------- loop over connectivity eigenvalues -------------------------
eigM = eig(M);        LLE = -Inf;
for mu = eigM.'
    RHS_poly = mu * hS * RHS_poly_base;
    pEq = polyadd(LHS_poly, -RHS_poly);
    lambdaRoots = roots(pEq);
    LLE = max(LLE, max(real(lambdaRoots)));
end

r0_analytic = r0;
LLE_analytic = LLE;

end

%% =====================================================================
%                       LOCAL FUNCTIONS
% =====================================================================

function s = polyadd(a,b)        % add two polynomial vectors
if numel(a)<numel(b), a=[zeros(1,numel(b)-numel(a)) a]; end
if numel(b)<numel(a), b=[zeros(1,numel(a)-numel(b)) b]; end
s = a+b;
end
%% ---------------------------------------------------------------------
%  Robust fixed‑point for *multi‑τ STD*
% ---------------------------------------------------------------------
function [r,b_mat,h,info] = fixed_point_multiSTD( ...
                u_ex,M,c_SFA,tau_a,F_STD,tau_b,tau_STD,n_a)

n   = numel(u_ex);      n_b = numel(tau_b);
gamma_mat = F_STD .* (tau_b(:)'/tau_STD);   %  n × n_b

% --- initial guess : tiny asymmetry to avoid clustered solution -------
r = max(0,u_ex) + 1e-3*randn(n,1);

% Anderson depth‑2 buffers
mAA  = 2;        Rbuf = zeros(n,mAA);  Gbuf = zeros(n,mAA);
tol_pic = 1e-4;  tol_final = 1e-10;
maxIt1  = 400;   maxIt2   = 200;

alpha = 1;  alpha_min = 1e-3;

F   = @(r) r - phi(r,u_ex,M,c_SFA,gamma_mat,tau_a,n_a);

% -------------- Stage‑1 : damped Picard + Anderson --------------------
for it = 1:maxIt1
    g   = F(r);     res = norm(g,inf);
    if res < tol_pic, break, end
    
    idx = mod(it-1,mAA)+1;
    Rbuf(:,idx) = r;   Gbuf(:,idx) = g;
    
    % Anderson extrapolation
    if it>=2
        dG = Gbuf(:,idx) - Gbuf(:,mod(idx-2,mAA)+1);
        dR = Rbuf(:,idx) - Rbuf(:,mod(idx-2,mAA)+1);
        beta = (dG(:)'*g(:)) / max(dG(:)'*dG(:),eps);
        rA   = r - beta*dR;
    else
        rA   = r - g;      % plain Picard
    end
    
    % monotone damping
    alpha = min(alpha*1.5,1);
    while true
        r_try = (1-alpha)*r + alpha*rA;
        if norm(F(r_try),inf) < res*(1-1e-4) || alpha<=alpha_min
            r = r_try;   break
        end
        alpha = alpha/2;
    end
end
it1 = it;  res1 = norm(F(r),inf);

% -------------- Stage‑2 : Newton–Krylov (GMRES) -----------------------
it2 = 0;
if res1 > tol_final
    for it2 = 1:maxIt2
        % auxiliaries
        [p,h,b_mat,dpdr,indA] = phb_from_r(r,gamma_mat,c_SFA,tau_a);
        g = F(r);  res = norm(g,inf);
        if res < tol_final, break, end
        
        % Jacobian‑vector product
        Afun = @(v) Jv_multiSTD(v,r,dpdr,M,c_SFA,n_a,indA);
        rhs  = -g;
        [dx,~,flag] = gmres(Afun,rhs,20,1e-8,20);
        if flag, warning('GMRES flag %d',flag); end
        
        % damped Newton line‑search
        alpha = 1;
        while true
            r_try = r + alpha*dx;
            if norm(F(r_try),inf) < res*(1-1e-4) || alpha<1e-6
                r = r_try; break
            end
            alpha = alpha/2;
        end
    end
end

% final auxiliaries
[p,h,b_mat] = phb_from_r(r,gamma_mat,c_SFA,tau_a);

% return info
info.it_picard = it1;
info.it_newton = it2;
info.it_total  = it1+it2;
info.res_final = norm(F(r),inf);
end
% ---------------------------------------------------------------------

% ---------- φ(r) = max(0,u_ex+M p - c_SFA Σa_k) ----------------------
function rF = phi(r,u_ex,M,c_SFA,gamma_mat,tau_a,n_a)
[p,~,~,~,~] = phb_from_r(r,gamma_mat,c_SFA,tau_a);
u_sfa = zeros(size(r));
for k = 1:n_a
    u_sfa = u_sfa + c_SFA .* r ./ (1 + c_SFA.*r/ (1/tau_a(k)) ); % Not quite right for multi-a SFA
end

% a_k = r / (1/tau_a_k + c_SFA_i*r). This is for single timescale.
% For multi-timescale SFA, the effective input is u_eff = u_d - sum(a_k)
% at fixed point a_k = c_SFA_k * r. so u_eff = u_d - (sum(c_SFA_k)) * r
% the c_SFA passed is a scalar strength, distributed over timescales.
% Original code used c_SFA .* (n_a*r)
u_eff = u_ex + M*p - c_SFA.*(n_a*r);

rF    = max(0,u_eff);
end

% ---------- p,h,b,dpdr from r  (multi‑τ) ------------------------------
function [p,h,b_mat,dpdr,indA] = phb_from_r(r,gamma_mat,c_SFA,tau_a)
[n,n_b] = size(gamma_mat);
n_a = numel(tau_a);
p   = zeros(n,1);
h   = zeros(n,1);
b_mat = zeros(n,n_b);
dpdr  = zeros(n,1);

for i = 1:n
    ri = r(i);
    gam = gamma_mat(i,:).';            % n_b×1
    if ri==0  % No STD/SFA if silent
        p(i)     = ri;
        if n_b>0, b_mat(i,:) = ones(1,n_b); end
        h(i)     = 1;
        dpdr(i)  = 1;
        continue
    end
    
    % --- Newton solve  f(p)=0 :  p - r ∏ 1/(1+γ p) -------------------
    pi = ri / (1+sum(gam)*ri);         % good initial guess
    for k = 1:50
        denom = 1 + gam*pi;
        if any(denom==0) % prevent division by zero
            pi = pi + 1e-9;
            denom = 1 + gam*pi;
        end
        hi    = prod(1./denom);        % h(p)
        f     = pi - ri*hi;
        if abs(f) < 1e-13, break, end
        hprime = hi * (-sum(gam./denom));
        dfdp   = 1 - ri*hprime;        %  = 1 + pi Σ γ/(1+γ p)
        pi     = pi - f/dfdp;
        pi     = max(pi,0);            % safeguard
    end
    p(i)   = pi;
    if n_b>0
        b_mat(i,:) = (1./(1+gam*pi)).';
        h(i)   = prod(b_mat(i,:));
        hprime = -h(i)*sum(gam./(1+gam*pi));
        dpdr(i)= h(i) / (1 - ri*hprime);   % chain rule
    else
        h(i) = 1;
        dpdr(i) = 1;
    end
end
indA = (r>0);
end

% % ---------- Jacobian‑vector product  J*v  (multi‑τ) --------------------
% function Jv = Jv_multiSTD(v,r,dpdr,M,c_SFA,n_a,indA)
% % F(r) = r - phi(r);    J = I - dphi/dr
% Jv = v;                               % start with identity
% w  = (v + c_SFA.*(n_a*v)) - M*(dpdr.*v); % Added parenthesis for clarity.
% ia = find(indA);                      % active neurons (ReLU slope 1)
% Jv(ia) = Jv(ia) - w(ia);
% end
% 
% % ---------- Jacobian‑vector product  J*v  (multi‑τ) --------------------
% function Jv = Jv_multiSTD(v,r,dpdr,M,c_SFA,n_a,indA)
% % F(r) = r - phi(r);    J = I - dphi/dr
% Jv = v;                               % start with identity
% w  = v + c_SFA.*(n_a*v) - M*(dpdr.*v);
% ia = find(indA);                      % active neurons (ReLU slope 1)
% Jv(ia) = w(ia);
% end

 function Jv = Jv_multiSTD(v,r,dpdr,M,c_SFA,n_a,indA)

    Jv = v;                               % identity part
    ia = find(indA);                      % active rows only
    Jv(ia) = v(ia) ...                    % I · v
             + c_SFA(ia).*n_a.*v(ia) ...  %  + c_SFA·n_a · v
             - M(ia,:)*(dpdr.*v);         %  – M·(dpdr ∘ v)
 end


