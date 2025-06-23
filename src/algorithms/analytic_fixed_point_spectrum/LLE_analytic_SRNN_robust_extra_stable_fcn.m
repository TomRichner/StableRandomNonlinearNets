function [r0_analytic, LLE_analytic] = LLE_analytic_SRNN_robust_extra_stable_fcn( ...
    n, nE, nI, M, u_ex_dc, ...
    n_a,  tau_a,  c_SFA, ...
    n_b,  tau_b,  F_STD, tau_STD, ...
    tau_d)
%% -----------------------------------------------------------
%  LLE_analytic_SRNN_robust_fcn
%  -----------------------------------------------------------
%  Analytic fixed-point + Lyapunov calculation for the SRNN rate
%  model with ReLU, multi-τ SFA and multi-τ STD.
%
%  Fixed point is obtained by **homotopy continuation** on DC:
%     σ : 0 → 1  in  9–20 sub-steps (adaptive bisection)
%     • damped Picard + Anderson (depth 2) to reach ‖F‖∞ < 1e-4
%     • Newton–Krylov (GMRES) to machine precision  (‖F‖∞ < 1e-11)
%
%  Lyapunov spectrum follows from the closed-form dispersion
%  relation evaluated on every eigen-value μ of the connectivity M.
%  -----------------------------------------------------------
%  Call signature identical to previous versions; just replace
%  the old .m file with this one.
%  -----------------------------------------------------------

%% 1 - make parameter vectors column-shaped
tau_a = tau_a(:);
tau_b = tau_b(:);

%% 2 - DC vector
u_ex_vec = u_ex_dc * ones(n,1);

%% 3 - fixed point (homotopy continuation version)
[r0,   b_mat, h0, infoFP] = fixed_point_multiSTD( ...
u_ex_vec, M, c_SFA, tau_a, F_STD, tau_b, tau_STD, n_a);

fprintf('FP | Picard %3d  Newton %3d   ‖g‖∞ = %.2e   max(r0)=%.4g\n', ...
infoFP.it_picard, infoFP.it_newton, infoFP.res_final, max(r0));

%% 4 - analytic Lyapunov spectrum  (use any active E-cell)
idx_active   = (r0>0);
iRef         = find(idx_active & (c_SFA>0 | F_STD>0), 1, 'first');
if isempty(iRef)
warning('No active excitatory neuron with adaptation – cannot compute analytic LLE.');
r0_analytic = r0;  LLE_analytic = NaN;  return
end

rS   = r0(iRef);
bS   = b_mat(iRef,:).';                 % n_b×1
hS   = prod(bS);
cS   = c_SFA(iRef);
Kvec = F_STD(iRef) * rS .* tau_b ./ (bS*tau_STD);

%% 4.1  build SFA polynomial  A(λ) = NumerA / DenA
DenA = 1;
for k = 1:n_a,  DenA = conv(DenA,[tau_a(k) 1]);  end
NumerA = DenA;
if cS>0
add = 0;
for k = 1:n_a
add = polyadd(add, cS*deconv(DenA,[tau_a(k) 1]));
end
NumerA = polyadd(NumerA, add);
end

%% 4.2  STD polynomial  B(λ) = NumerB / DenB
DenB = 1;
for m = 1:n_b,  DenB = conv(DenB,[tau_b(m) 1]);  end
NumerB = DenB;
for m = 1:n_b
NumerB = polyadd(NumerB, Kvec(m)*deconv(DenB,[tau_b(m) 1]));
end

%% 4.3  dendrite  C(λ)=τ_d λ+1   and dispersion polynomial
Cpoly    = [tau_d 1];
LHS_poly = conv(conv(Cpoly,NumerA),DenB);     % C · NumerA · DenB
RHS_base = conv(NumerB,DenA);                 % NumerB · DenA

eigM  = eig(M);
LLE   = -Inf;
for mu = eigM.'
pEq = polyadd(LHS_poly, -mu*hS*RHS_base);
LLE = max(LLE, max(real(roots(pEq))));
end

%% 5 - outputs
r0_analytic = r0;
LLE_analytic = LLE;
end  % ================= END MAIN FUNCTION ==============================


%% =====================================================================
%  Helper: add two polynomials (highest power first)
% =====================================================================
function s = polyadd(a,b)
if numel(a)<numel(b), a=[zeros(1,numel(b)-numel(a)) a]; end
if numel(b)<numel(a), b=[zeros(1,numel(a)-numel(b)) b]; end
s = a+b;
end

%% =====================================================================
%  Fixed-point solver with DC-continuation
% =====================================================================
function [r,b_mat,h,info] = fixed_point_multiSTD( ...
  u_ex_target, M, c_SFA, tau_a, ...
  F_STD,      tau_b, tau_STD, n_a)

n    = numel(u_ex_target);
n_b  = numel(tau_b);
gamma= F_STD .* (tau_b(:)'/tau_STD);           % n×n_b

% --- continuation ladder in σ ------------------------------------------------
sigma = [0 .02 .05 .1 .2 .4 .6 .8 1];          % will adaptively bisect

% start from silent root
r = zeros(n,1);

for s = 1:numel(sigma)
u_ex = sigma(s)*u_ex_target;

[r,b_mat,h,itP,itN,res] = solve_single_level(r,u_ex,M,c_SFA, ...
                                     gamma,tau_a,n_a);
% adaptive refinement
if res>1e-4 && s<numel(sigma)
sigma = sort([sigma, (sigma(s)+sigma(s+1))/2]);
continue
end
end

info.it_picard = itP;
info.it_newton = itN;
info.res_final = res;
end

% ---------- solve at a single continuation step -----------------------
function [r,b_mat,h,it_pic,it_new,res_fin] = solve_single_level( ...
       r0,u_ex,M,c_SFA,gamma_mat,tau_a,n_a)

n        = numel(u_ex);
tol_pic  = 1e-4;     tol_fin = 1e-11;
maxItPic = 400;      maxItNew = n;

F  = @(x) x - phi(x,u_ex,M,c_SFA,gamma_mat,n_a);

% ---------- Picard + Anderson depth-2 ---------------------------------
r     = r0;     mAA=2;  Rbuf=zeros(n,mAA); Gbuf=Rbuf;
for it_pic = 1:maxItPic
g   = F(r);   res = norm(g,inf);
if res<tol_pic, break, end
idx = mod(it_pic-1,mAA)+1;
Rbuf(:,idx)=r;   Gbuf(:,idx)=g;
if it_pic>=2
dG = Gbuf(:,idx)-Gbuf(:,mod(idx-2,mAA)+1);
dR = Rbuf(:,idx)-Rbuf(:,mod(idx-2,mAA)+1);
beta = (dG(:)'*g(:))/max(dG(:)'*dG(:),eps);
rA   = r - beta*dR;
else
rA   = r - g;
end
alpha=1;
while true
r_try=(1-alpha)*r+alpha*rA;
if norm(F(r_try),inf)<res*(1-1e-4) || alpha<1e-3
r=r_try; break
end
alpha=alpha/2;
end
end

% ---------- Newton–Krylov ---------------------------------------------
for it_new = 0:maxItNew
[p,h,b_mat,dpdr,indA] = phb_from_r(r,gamma_mat);
g     = F(r);     res_fin = norm(g,inf);
if res_fin<tol_fin, break, end

Af   = @(v) Jv_multi(v,dpdr,M,c_SFA,n_a,indA);
rhs  = -g;
% restart = min(20,n);
restart = n;
[dx,flag] = gmres(Af,rhs,restart,1e-8,n);
if flag, warning('GMRES flag=%d',flag); end
alpha=1;
while true
r_try=r+alpha*dx;
if norm(F(r_try),inf)<res_fin*(1-1e-4) || alpha<1e-6
r=r_try; break
end
alpha=alpha/2;
end
end
end  % solve_single_level

%% ---------- φ(r) ------------------------------------------------------
function rF = phi(r,u_ex,M,c_SFA,gamma_mat,n_a)
[p,~,~,~,~] = phb_from_r(r,gamma_mat);          % STD
u_eff = u_ex + M*p - c_SFA.*(n_a*r);            % SFA feedback
rF    = max(0,u_eff);
end

%% ---------- p,h,b,dpdr  ----------------------------------------------
function [p,h,b_mat,dpdr,indA] = phb_from_r(r,gamma_mat)
[n,n_b]=size(gamma_mat);
p=zeros(n,1); h=ones(n,1); b_mat=ones(n,n_b); dpdr=ones(n,1);
for i=1:n
ri=r(i); gam=gamma_mat(i,:)';
if ri==0 || all(gam==0), continue, end
pi = ri/(1+ri*sum(gam));                 % init
for k=1:30
denom = 1+gam*pi;  hi=prod(1./denom);
f = pi - ri*hi;   if abs(f)<1e-13, break, end
dfdp = 1 - ri*hi*(-sum(gam./denom));
pi = max(pi - f/dfdp, 0);
end
p(i)=pi; b_mat(i,:)=1./(1+gam*pi)'; h(i)=prod(b_mat(i,:));
hprime = -h(i)*sum(gam./(1+gam*pi));
dpdr(i)= h(i)/(1 - ri*hprime);
end
indA = r>0;
end

%% ---------- Jacobian-vector product ----------------------------------
function Jv = Jv_multi(v,dpdr,M,c_SFA,n_a,indA)
Jv = v;
ia = find(indA);
Jv(ia) = v(ia) + c_SFA(ia).*n_a.*v(ia) - M(ia,:)*(dpdr.*v);
end
