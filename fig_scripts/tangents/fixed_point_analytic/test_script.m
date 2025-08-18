% SRNN_fp_lle_compare.m
close all; clc; clear SRNN_NL; clearvars -except seed;

%% Network & params (mirrors SRNN_example_tseries)
seed = 107; rng(seed,'twister');
n = 10;
EI_frac = 0.7;
mean_in_out_degree = 5;
density  = mean_in_out_degree/(n-1);
sparsity = 1 - density;

% Weight scales (like your example)
scale = 0.5/0.79782; 
w.EE=scale*1.00; 
w.EI=scale*1; 
w.IE=scale*1; 
w.II=scale*0.5;
w.selfE = 0;
w.selfI = 0;
[M, EI_vec] = generate_M_no_iso(n, w, sparsity, EI_frac);
isE = EI_vec(:)>0; isI = ~isE;

% Time constants
tau_d   = 0.025;    % s
tau_STD = 1.0;    % s
K       = 3;      % SFA timescales (E only)
tau_a_E = logspace(log10(0.3), log10(15), K).'; % s
tau_b   = tau_STD;   % s

% SFA per-neuron coeff vector (matches your example convention)
c_SFA = zeros(n,1); c_SFA(isE) = 0.5/K;  % so total SFA gain per E neuron is 0.5
% STD gating (E only)
F_STD = zeros(n,1); F_STD(isE) = 1;

% Package params for SRNN_NL
E_idx = find(isE); I_idx = find(isI);
params = package_params(nnz(isE), nnz(isI), E_idx, I_idx, ...
                        K, 0, 1, 0, ...
                        tau_a_E(:).', [], tau_b, [], ...
                        tau_d, n, M, c_SFA, F_STD, tau_STD, EI_vec);
params.activation_function = @(x) max(x,0); % ReLU

%% Fixed-point
DC = 0.1;  % constant input to all neurons
opts_fp = struct('tol',1e-12,'maxit',5e4,'alpha',0.6,'Fa',ones(K,1),'verbose',true);
[u_fp, A_fp, b_fp, r_fp, info_fp] = SRNN_fixed_point( ...
    M, isE, c_SFA, tau_a_E(:), tau_b, tau_d, tau_STD, DC, F_STD, opts_fp);

% Sanity: algebraic residuals
res_u = norm(u_fp - (DC + M*(r_fp.*b_fp)), inf);
res_a = norm(A_fp(isE,:) - r_fp(isE).*ones(1,K), inf);  % should be 0 if Fa=1
res_b = norm(b_fp(isE) - tau_STD./(tau_STD + tau_b*max(r_fp(isE),0)), inf);
fprintf('\nAlgebraic residuals:  ||u-DC-M(rb)||_∞=%.3g,  ||a-r||_∞=%.3g,  ||b-b*(r)||_∞=%.3g\n',res_u,res_a,res_b);

%% Analytic LLE at FP
[J_fp, lambda_max_fp, gamma_fp] = SRNN_jacobian_from_state( ...
    M, isE, c_SFA, tau_a_E(:), tau_b, tau_d, tau_STD, ...
    u_fp, A_fp, b_fp, r_fp, F_STD);
fprintf('Analytic LLE at fixed point:   %.6f  (1/s)\n', lambda_max_fp);

%% Simulation with constant DC
fs = 1000; dt = 1/fs;  T = [-10 80];
nt = round((T(2)-T(1))*fs)+1; t = linspace(T(1),T(2),nt).';
u_ex = DC*ones(n,nt);

a0_E = zeros(nnz(isE)*K,1); b0_E = ones(nnz(isE)*1,1);
X0 = [a0_E; zeros(0,1); b0_E; zeros(0,1); zeros(n,1)];  % [a_E; a_I; b_E; b_I; u_d]
SRNN_wrapper = @(tt,XX) SRNN_NL(tt,XX,t,u_ex,params);
ode_options = odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',dt,'InitialStep',0.1*dt);

fprintf('--- Running simulation (ode15s) ---\n');
[t_out, X] = ode15s(SRNN_wrapper, t, X0, ode_options);
assert(all(abs(t_out - t) < 1e-12));

[a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts] = unpack_SRNN_state(X, params);
[r_ts, p_ts] = compute_dependent_variables(a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params);

% Late-time averages
tail = t >= (T(2)-1.0);
u_mean = mean(u_d_ts(:,tail),2);
r_mean = mean(r_ts(:,tail),2);
fprintf('\n||u_mean - u_fp||_∞ = %.3g\n', norm(u_mean - u_fp, inf));
fprintf('||r_mean - r_fp||_∞ = %.3g\n', norm(r_mean - r_fp, inf));

r_comp = [r_ts(:,end), r_fp]

%% Benettin LLE
T_lya_1 = -5; lya_dt = 0.5*tau_d; d0 = 1e-3;
idx0 = find(t>=T_lya_1,1,'first');
[LLE_ben, local_lya, finite_lya, t_lya] = benettin_algorithm( ...
    X(idx0:end,:), t(idx0:end), dt, fs, d0, T, lya_dt, params, ode_options, ...
    @SRNN_NL, t, u_ex, @ode15s);
fprintf('\nBenettin LLE (trajectory):     %.6f  (1/s)\n', LLE_ben);

%% Jacobian LLE at final simulated state
u_end = u_d_ts(:,end);
A_end = zeros(n,K);  if K>0, A_end(isE,:) = a_E_ts(:,:,end); end
b_end = ones(n,1);   b_end(isE) = squeeze(b_E_ts(:,:,end));  % n_b_E=1 assumed
r_end = r_ts(:,end);

[J_end, lambda_max_end, gamma_end] = SRNN_jacobian_from_state( ...
    M, isE, c_SFA, tau_a_E(:), tau_b, tau_d, tau_STD, ...
    u_end, A_end, b_end, r_end, F_STD);
fprintf('Analytic LLE at final state:   %.6f  (1/s)\n', lambda_max_end);

%% Plot finite-time LLE vs analytic values
figure; plot(t_lya, local_lya, 'LineWidth',1); hold on;
yline(lambda_max_fp,'k--','FP LLE'); 
yline(lambda_max_end,'r-.','Final-state LLE');
xlabel('time (s)'); ylabel('finite-time LLE (1/s)');
title('Benettin finite-time LLE vs analytic'); grid on;
legend('finite-time LLE','FP LLE','final-state LLE');


%% Plot firing rates
figure;
hold on;
hE = plot(t, r_ts(isE,:).');
hI = plot(t, r_ts(isI,:).', 'r');
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Simulated Firing Rates');
if ~isempty(hE) && ~isempty(hI)
    legend([hE(1), hI(1)], 'Excitatory', 'Inhibitory');
elseif ~isempty(hE)
    legend(hE(1), 'Excitatory');
elseif ~isempty(hI)
    legend(hI(1), 'Inhibitory');
end
grid on;


