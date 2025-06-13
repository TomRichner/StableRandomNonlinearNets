%% -----------------------------------------------------------
%  File:   LLE_analytic_SRNN_stable.m
%  Author: ChatGPT (OpenAI o3)
%  Date:   2025‑06‑12
%
%  Fixed‑point search is now robust against max‑flip oscillations
%  using damped updates + Anderson acceleration (depth = 2).
%  -----------------------------------------------------------

clear;  clc;

%% ---------------- 1.  USER PARAMETERS  (unchanged) --------------------

n     = 10;      fracE = 0.7;   nE = round(fracE*n);  nI = n-nE;

rng(7);

g = 0.5;  
density = 0.5;  
M = g/sqrt(n) * (randn(n) .* (rand(n)<density));
M(:,1:nE)      =  abs(M(:,1:nE));
M(:,nE+1:end)  = -abs(M(:,nE+1:end));
M(eye(size(M),'logical')) = 0; % main diag is zero

u_ex = 0.1*ones(n,1);

n_a   = 3;     tau_a = logspace(log10(0.3),log10(15),n_a)';  % column
c_SFA = ones(n,1);  c_SFA(nE+1:end)=0;

tau_b = 4;     tau_STD = 0.5;
F_STD = ones(n,1);  F_STD(nE+1:end)=0;

tau_d = 0.025;

%% ---------------- 2.  FIXED‑POINT via damped Picard+AA ---------------

[r0,b0,h0,itFP] = fixed_point_SRNN(u_ex,M,c_SFA,tau_a,...
                                   F_STD,tau_b,tau_STD,n_a);
fprintf('FP converged in %d iterations   max(r0)=%.4g Hz\n',itFP,max(r0));

%% ---------------- 3.  ANALYTIC LLE  (exactly as before) --------------

idxA = (r0>0);                   % active neurons for analytic block
if any(diff(c_SFA(idxA))) || any(diff(F_STD(idxA)))
    warning('heterogeneous adaptation parameters – spectrum approximation used');
end
% Pick *representative* active E–cell parameters (all identical here)
iRef   = find(idxA & (c_SFA>0),1,'first');
cS     = c_SFA(iRef);      rS = r0(iRef);
bS     = b0(iRef);         hS = h0(iRef);
gammaS = F_STD(iRef)*(tau_b/tau_STD);

% SFA polynomials
DenA = 1;
for k=1:n_a,  DenA = conv(DenA,[tau_a(k) 1]);  end
NumerA = DenA;
for k=1:n_a
    NumerA = polyadd(NumerA, cS*deconv(DenA,[tau_a(k) 1]));
end

% STD polynomial  (one τ_b)
K      = F_STD(iRef)*rS/(bS*tau_STD);
DenB   = [1 1/tau_b];      NumerB = [1 1/tau_b+K];

Cpoly  = [tau_d 1];        % λ τ_d + 1
ABCnum = conv(conv(Cpoly,NumerA),NumerB);
ABDden = conv(DenA,DenB);

eigM   = eig(M);           LLE = -Inf;
for mu = eigM.'
    pEq = polyadd( ABCnum , -mu*hS*ABDden );
    rootsMu = roots(pEq);
    LLE = max(LLE, max(real(rootsMu)));
end
fprintf('Analytic λ_max = %+8.5f  1/s\n',LLE);

%% ---------------- 4.  helper:  polyadd -------------------------------
function s = polyadd(a,b)
if numel(a)<numel(b), a=[zeros(1,numel(b)-numel(a)) a]; end
if numel(b)<numel(a), b=[zeros(1,numel(a)-numel(b)) b]; end
s = a+b;
end

%% =====================================================================
%  SUBFUNCTION :  robust fixed‑point finder
%  =====================================================================
function [r,b,h,it] = fixed_point_SRNN(u_ex,M,c_SFA,tau_a,...
                                       F_STD,tau_b,tau_STD,n_a)

n   = numel(u_ex);
gamma = F_STD.*(tau_b/tau_STD);

% Anderson‑acceleration buffers (depth = 2)
depth = 2;  R = zeros(n,depth);  G = zeros(n,depth);

r  = max(0,u_ex);                % initial guess
tol = 1e-6;   maxIt = 10000;
eta = 1;        eta_min = 1e-3;  eta_max = 1;

for it = 1:maxIt
    
    % ---- forward map F(r) --------------------------------------------
    p = r./(1+gamma.*r);
    u_eff = u_ex + M*p - c_SFA.*(n_a*r);
    rF = max(0,u_eff);           % F(r)
    
    % residual
    g  = rF - r;
    if norm(g,inf) < tol,  break;  end
    
    % ---------- Anderson mix (depth 2) --------------------------------
    idx = mod(it-1,depth)+1;
    R(:,idx) = r;  G(:,idx) = g;
    m = min(it,depth);
    
    if m>1
        dG = G(:,idx) - G(:,mod(idx-2,depth)+1);
        dR = R(:,idx) - R(:,mod(idx-2,depth)+1);
        beta = (dG(:)'*g(:)) / (dG(:)'*dG(:)+eps);
        rA   = r - beta*dR;
    else
        rA   = r + g;            % plain Picard first step
    end
    
    % ---------- damping / back‑tracking -------------------------------
    eta = min(eta*1.2, eta_max);      % try to increase cautiously
    while true
        r_new = r + eta*(rA - r);     % candidate
        p_new = r_new./(1+gamma.*r_new);
        u_eff_new = u_ex + M*p_new - c_SFA.*(n_a*r_new);
        rF_new = max(0,u_eff_new);
        g_new  = rF_new - r_new;
        if norm(g_new,inf) < (1-1e-4)*norm(g,inf) || eta<=eta_min
            r = r_new;   break
        end
        eta = eta/2;                % back‑track
    end
end

if it==maxIt,  warning('fixed‑point solver reached maxIt'); end

% post‑solve: STD variables
p = r./(1+gamma.*r);    h = 1./(1+gamma.*r);  b = h;
end
