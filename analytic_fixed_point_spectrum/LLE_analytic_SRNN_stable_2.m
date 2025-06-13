%% -----------------------------------------------------------
%  File:   LLE_analytic_SRNN_robust.m
%  Author: ChatGPT (OpenAI o3)  –  2025‑06‑12
%
%  Fixed‑point search:
%     • Stage‑1: damped Picard + Anderson depth‑2
%     • Stage‑2: Newton–Krylov (GMRES)
%
%  Works for n up to a few‑thousand on a laptop.
%  -----------------------------------------------------------

clear;  clc;

%% ---------------- 1.  USER PARAMETERS  -------------------------------

n       = 10;                  % try a few hundred neurons
fracE   = 0.7;   nE = round(fracE*n);  nI = n-nE;

g       = 0.8;    density = 0.55;
% rng(7);
M = g/sqrt(n) .* (randn(n) .* (rand(n)<density));
M(:,1:nE)      =  abs(M(:,1:nE));   % E columns +
M(:,nE+1:end)  = -abs(M(:,nE+1:end)); % I columns –
M(eye(size(M),'logical')) = 0; % main diag is zero

u_ex   = 0.05*ones(n,1);            % constant DC

% adaptation
n_a   = 3;  tau_a = logspace(log10(0.3),log10(6),n_a)';
c_SFA = ones(n,1);  c_SFA(nE+1:end)=0;

tau_b = 4;   tau_STD = 0.5;
F_STD = ones(n,1);  F_STD(nE+1:end)=0;

tau_d = 0.025;

%% ---------------- 2.  FIXED‑POINT  -----------------------------------

[r0,b0,h0,infoFP] = fixed_point_SRNN(u_ex,M,c_SFA,tau_a,...
                                     F_STD,tau_b,tau_STD,n_a);

fprintf('FP  |  it=%d  ‖g‖∞=%.2e   stage‑2 it=%d   max(r0)=%.4g Hz\n',...
         infoFP.it_total,infoFP.res_final,infoFP.it_newton,max(r0));

%% ---------------- 3.  ANALYTIC LLE  ----------------------------------

idxA  = (r0>0);
iRefE = find(idxA & (c_SFA>0),1,'first');   % first active E neuron

cS  = c_SFA(iRefE);     rS = r0(iRefE);
bS  = b0(iRefE);        hS = h0(iRefE);
gammaS = F_STD(iRefE)*(tau_b/tau_STD);

% --- build scalar polynomials (as in the earlier answer) ---------------
DenA = 1;
for k=1:n_a,  DenA = conv(DenA,[tau_a(k) 1]); end
NumerA = DenA;
for k=1:n_a
    NumerA = polyadd(NumerA, cS*deconv(DenA,[tau_a(k) 1]));
end

K      = F_STD(iRefE)*rS/(bS*tau_STD);
DenB   = [1 1/tau_b];      NumerB = [1 1/tau_b+K];

Cpoly  = [tau_d 1];
ABCnum = conv(conv(Cpoly,NumerA),NumerB);
ABDden = conv(DenA,DenB);

eigM   = eig(M);           LLE = -Inf;
for mu = eigM.'
    pEq = polyadd(ABCnum, -mu*hS*ABDden);
    lamb = roots(pEq);
    LLE  = max(LLE, max(real(lamb)));
end
fprintf('Analytic Λ_max = %+8.5f  1/s\n',LLE);

%% ---------------- helper: polyadd ------------------------------------
function s = polyadd(a,b)
if numel(a)<numel(b), a=[zeros(1,numel(b)-numel(a)) a]; end
if numel(b)<numel(a), b=[zeros(1,numel(a)-numel(b)) b]; end
s = a+b;
end

%% =====================================================================
%  SUBFUNCTION :  Robust fixed‑point finder
%  =====================================================================
function [r,b,h,info] = fixed_point_SRNN(u_ex,M,c_SFA,tau_a,...
                                         F_STD,tau_b,tau_STD,n_a)
n = numel(u_ex);
gamma = F_STD .* (tau_b/tau_STD);

% ---------- initial guess: tiny noise to break symmetry ---------------
r = max(0, u_ex + 1e-3*randn(n,1));

% Anderson depth
mAA = 2;
Rbuf = zeros(n,mAA);  Gbuf = zeros(n,mAA);

tol_pic  = 1e-4;      % switch to Newton when ‖g‖∞ < tol_pic
tol_final = 1e-10;    % final required accuracy
maxIt1 = 300;         % Picard/AA iterations
maxIt2 = 30;          % Newton iterations

alpha = 1;  alpha_min = 1e-3;

% residual function handle ----------------------------------------------
F = @(r) r - phi(r,u_ex,M,c_SFA,gamma,n_a);

% ---------- Stage‑1 : damped Picard + Anderson -------------------------
for it = 1:maxIt1
    g  = F(r);
    res = norm(g,inf);
    if res < tol_pic, break, end
    
    idx = mod(it-1,mAA)+1;
    Rbuf(:,idx) = r;
    Gbuf(:,idx) = g;
    
    % Anderson step ------------------------------------------------------
    if it>=2
        dG = Gbuf(:,idx) - Gbuf(:,mod(idx-2,mAA)+1);
        dR = Rbuf(:,idx) - Rbuf(:,mod(idx-2,mAA)+1);
        beta = (dG(:)'*g(:)) / max(dG(:)'*dG(:),eps);
        rA = r - beta*dR;
    else
        rA = r - g;    % simple Picard on first step
    end
    
    % monotone back‑tracking on residual --------------------------------
    alpha = min(alpha*1.5, 1);
    while true
        r_try = (1-alpha)*r + alpha*rA;
        res_try = norm(F(r_try),inf);
        if res_try < res*(1-1e-4) || alpha<=alpha_min
            r = r_try;   break
        end
        alpha = alpha/2;
    end
end
it1 = it;   res1 = norm(F(r),inf);

% ---------- Stage‑2 : Newton–Krylov if needed --------------------------
it2 = 0;
if res1 > tol_final
    Jv = @(v,ri,pi,hi,indA) jacobian_vec(v,ri,pi,hi,indA,M,gamma,c_SFA,n_a);
    for it2 = 1:maxIt2
        % recompute auxiliaries
        [p,h,indA] = ph_from_r(r,gamma);
        g = F(r);    res = norm(g,inf);
        if res < tol_final, break, end
        
        % GMRES(20) to solve J * dx = -g
        rhs  = -g;
        Afun = @(v) Jv(v,r,p,h,indA);   % J*v
        [dx,~,gmflag] = gmres(Afun,rhs,20,1e-8,20);
        if gmflag, warning('GMRES did not fully converge'); end
        
        % damped Newton line‑search (Armijo)
        alpha = 1;
        while true
            r_try = r + alpha*dx;
            res_try = norm(F(r_try),inf);
            if res_try < res*(1-1e-4) || alpha <= 1e-6
                r = r_try;  break
            end
            alpha = alpha/2;
        end
    end
end

% ----------- final auxiliaries -----------------------------------------
[p,h,~] = ph_from_r(r,gamma);   b = h;

info.it_total  = it1+it2;
info.it_picard = it1;
info.it_newton = it2;
info.res_final = norm(F(r),inf);
end

% ---------- φ : ReLU map -----------------------------------------------
function rF = phi(r,u_ex,M,c_SFA,gamma,n_a)
[p,~,~] = ph_from_r(r,gamma);
u_eff = u_ex + M*p - c_SFA.*(n_a*r);
rF = max(0,u_eff);
end

% ---------- compute p,h and active‑set flag ----------------------------
function [p,h,indA] = ph_from_r(r,gamma)
den = 1 + gamma.*r;
p = r ./ den;
h = 1 ./ den;
indA = r>0;
end

% ---------- Jacobian‑vector product  J*v  (no explicit J) --------------
function Jv = jacobian_vec(v,r,p,h,indA,M,gamma,c_SFA,n_a)
% For silent neurons derivative of ReLU is zero → row of identity.
Jv = v;                       % start with I*v
% active rows only:
ia = find(indA);
if isempty(ia), return, end

% components that enter through ReLU argument
w   = (1 + c_SFA.*n_a).*v;    % diagonal part

% term −M * ( d p / d r  * v )
dpdr = 1./(1+gamma.*r).^2;    % element‑wise derivative
w = w - M*(dpdr.*v);

Jv(ia) = w(ia);
end
