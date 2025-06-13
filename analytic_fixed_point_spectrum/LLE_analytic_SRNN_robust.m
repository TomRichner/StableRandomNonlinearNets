%% -----------------------------------------------------------
%  File :  LLE_analytic_SRNN_multiSTD.m
%  Author: ChatGPT (OpenAI o3)
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

clear;  clc;

%% ---------------- 1.  USER‑EDITABLE PARAMETERS ------------------------

n        = 10;                % total neurons
fracE    = 0.8;  nE = round(fracE*n);  nI = n-nE;

% Connectivity -----------------------------------------------------------
g        = 0.8;                % weight scale
density  = 0.55;               % connection probability
rng(42);
M = g/sqrt(n) .* (randn(n) .* (rand(n)<density));
M(:,1:nE)      =  abs(M(:,1:nE));      % E→· positive
M(:,nE+1:end)  = -abs(M(:,nE+1:end));  % I→· negative (Dale)
M(eye(size(M),'logical')) = 0; % main diag is zero

% External drive (constant DC) ------------------------------------------
u_ex     = 0.2*ones(n,1);

% SFA  ---------------------------------------------------------------
n_a      = 3;                                % # SFA timescales
tau_a    = logspace(log10(0.3),log10(15),n_a)';  % column (s)
c_SFA    = ones(n,1);           c_SFA(nE+1:end)=0;  % none on I

% *** STD :   now allow >1 time‑constant *****************************
n_b = 2;
tau_b    = logspace(log10(0.6),log10(9),n_b);            %  <-- EDIT HERE (s)   length = n_b

tau_STD  = 0.5;                 % release‑probability recovery (s)
F_STD    = ones(n,1);           F_STD(nE+1:end)=0;   % none on I

% Dendritic low‑pass -----------------------------------------------------
tau_d    = 0.025;               % (s)

%% ---------------- 2.  FIXED‑POINT SEARCH ------------------------------

[r0,b0_mat,h0,infoFP] = fixed_point_multiSTD( ...
        u_ex,M,c_SFA,tau_a,F_STD,tau_b,tau_STD,n_a);

fprintf('FP  |  Picard it=%d  Newton it=%d  ‖g‖∞=%.1e   max(r0)=%.4g\n',...
        infoFP.it_picard,infoFP.it_newton,infoFP.res_final,max(r0));

%% ---------------- 3.  ANALYTIC LLE  (multi‑τ STD) --------------------

% Choose a representative *active excitatory* neuron
idxA = (r0>0);
iRef = find(idxA & (c_SFA>0),1,'first');
if isempty(iRef)
    error('No active excitatory neuron – increase DC or weights.');
end

% Shorthand for that neuron
rS  = r0(iRef);
bS  = b0_mat(iRef,:).';                 % n_b × 1
hS  = h0(iRef);                         % product of bS
cS  = c_SFA(iRef);
gamma_vec = F_STD(iRef)*tau_b(:)/tau_STD;   % n_b ×1
Kvec = F_STD(iRef) * rS ./ (bS * tau_STD);  % K_m  (n_b×1)

% ---------- polynomials :   SFA (A)  and  STD (B) ----------------------
DenA = 1;
for k = 1:n_a
    DenA = conv(DenA,[tau_a(k) 1]);            % ∏ (τ_a λ+1)
end
NumerA = DenA;
for k = 1:n_a
    NumerA = polyadd(NumerA, cS*deconv(DenA,[tau_a(k) 1]));
end

% multi‑τ  STD
DenB = 1;
for m = 1:n_b
    DenB = conv(DenB,[1 1/tau_b(m)]);          % ∏ (λ+1/τ_bm)
end
NumerB = DenB;
for m = 1:n_b
    Partial = deconv(DenB,[1 1/tau_b(m)]);
    NumerB  = polyadd(NumerB, Kvec(m)*Partial);
end

% dendrite factor  (λ τ_d + 1)
Cpoly  = [tau_d 1];
ABCnum = conv(conv(Cpoly,NumerA),NumerB);
ABDden = conv(DenA,DenB);

% ---------- loop over connectivity eigenvalues -------------------------
eigM = eig(M);        LLE = -Inf;
for mu = eigM.'
    pEq = polyadd(ABCnum, -mu*hS*ABDden);     % dispersion poly
    lambdaRoots = roots(pEq);
    LLE = max(LLE, max(real(lambdaRoots)));
end
fprintf('Analytic  Λ_max = %+8.5f  1/s\n',LLE);

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
maxIt1  = 400;   maxIt2   = 40;

alpha = 1;  alpha_min = 1e-3;

F   = @(r) r - phi(r,u_ex,M,c_SFA,gamma_mat,n_a);

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
        [p,h,b_mat,dpdr,indA] = phb_from_r(r,gamma_mat);
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
[p,h,b_mat] = phb_from_r(r,gamma_mat);

% return info
info.it_picard = it1;
info.it_newton = it2;
info.it_total  = it1+it2;
info.res_final = norm(F(r),inf);
end
% ---------------------------------------------------------------------

% ---------- φ(r) = max(0,u_ex+M p - c_SFA n_a r) ----------------------
function rF = phi(r,u_ex,M,c_SFA,gamma_mat,n_a)
[p,~,~,~,~] = phb_from_r(r,gamma_mat);
u_eff = u_ex + M*p - c_SFA.*(n_a*r);
rF    = max(0,u_eff);
end

% ---------- p,h,b,dpdr from r  (multi‑τ) ------------------------------
function [p,h,b_mat,dpdr,indA] = phb_from_r(r,gamma_mat)
[n,n_b] = size(gamma_mat);
p   = zeros(n,1);
h   = zeros(n,1);
b_mat = zeros(n,n_b);
dpdr  = zeros(n,1);

for i = 1:n
    ri = r(i);
    gam = gamma_mat(i,:).';            % n_b×1
    if all(gam==0) || ri==0            % inhibitory or silent
        p(i)     = ri;
        b_mat(i,:) = ones(1,n_b);
        h(i)     = 1;
        dpdr(i)  = 1;
        continue
    end
    
    % --- Newton solve  f(p)=0 :  p - r ∏ 1/(1+γ p) -------------------
    pi = ri / (1+sum(gam)*ri);         % good initial guess
    for k = 1:50
        denom = 1 + gam*pi;
        hi    = prod(1./denom);        % h(p)
        f     = pi - ri*hi;
        if abs(f) < 1e-13, break, end
        hprime = hi * (-sum(gam./denom));
        dfdp   = 1 - ri*hprime;        %  = 1 + pi Σ γ/(1+γ p)
        pi     = pi - f/dfdp;
        pi     = max(pi,0);            % safeguard
    end
    p(i)   = pi;
    b_mat(i,:) = (1./(1+gam*pi)).';
    h(i)   = prod(b_mat(i,:));
    
    hprime = -h(i)*sum(gam./(1+gam*pi));
    dpdr(i)= h(i) / (1 - ri*hprime);   % chain rule
end
indA = (r>0);
end

% ---------- Jacobian‑vector product  J*v  (multi‑τ) --------------------
function Jv = Jv_multiSTD(v,r,dpdr,M,c_SFA,n_a,indA)
% F(r) = r - phi(r);    J = I - dphi/dr
Jv = v;                               % start with identity
w  = v + c_SFA.*(n_a*v) - M*(dpdr.*v);
ia = find(indA);                      % active neurons (ReLU slope 1)
Jv(ia) = w(ia);
end
