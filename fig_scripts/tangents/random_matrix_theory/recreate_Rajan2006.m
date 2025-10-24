close all
clear
clc

set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 18);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

n = 500;
b = 1;
mu = 0;

% mean_in_out_degree = 500; % desired mean number of connections in and out
mean_in_out_degree = 50; % desired mean number of connections in and out
density = mean_in_out_degree/(n) % each neuron can make up to n-1 connections with other neurons

% disk_radius=b*sqrt(density)
disk_radius = sqrt(n)*b*sqrt(density)
% disk_radius = sqrt(n)*b;
% % disk_radius = b*sqrt(density);

sigma = disk_radius;
% sigma = 1;

%% Generate A, Ginibre

A = b*randn(n,n);
% A = b*randn(n,n)./sqrt(n) + mu;

% A = b*randn(n,n) + mu;


%% weakly enforce Dale's like Rajan 2006
f = 0.5; % fraction excitatory
isE = zeros(n,1,'logical'); % is excitatory logical vector
isE(1:round(n*f),1) = true;

mu_E = 0.5; % mean of excitatory 
mu_I = -f*mu_E./(1-f) % mean of inhibitory, depends on mu_E to ensure f*mu_E+(1-f)*mu_I = 0 for "balance"

A(:, isE) =  A(:, isE) + f*mu_E;   % shifts excitatory mean positive
A(:, not(isE)) = A(:, not(isE)) + f*mu_I;   % shifts inhibitory mean negative

%% make it sparse
sparse_mask = rand(n,n)<(1-density);
A(sparse_mask) = 0; % make it sparse
% A(eye(size(A),'logical')) = 0; % no self connections, especially if Dales

%% check A stats

actual_var_A = var(A(:))
actual_stdev_A = std(A(:))
actual_mean_A = mean(A(:))
actual_mean_in_out_degree = sum(abs(A(:))>0)/n
actual_density = sum(abs(A(:))>0)/n^2

%% enforce row sums

% % fix only half of the rows -- doesn't helf
% for i_ro = 1:fix(0.5*n) 
% 
%     A(i_ro,:) = A(i_ro,:) - sum(A(i_ro,:))./n;
% 
% end

% fix all rows
for i_ro = 1:n % fix all

    A(i_ro,:) = A(i_ro,:) - sum(A(i_ro,:))./n;

end


%% Correct A onto edge of stability by shifting
corrected_A = A-sigma*eye(size(A));
corrected_var_A = var(A(:))
corrected_stdev_A = std(A(:))
corrected_mean_A = mean(A(:))

eig_A = eig(corrected_A);


%% make plots

% plot eigenvalues
figure(1)
set(gcf,'Position',[1965         897         560         420])
t = 0:.01:2*pi;
% plot(sin(t)-sigma,cos(t),'g','LineWidth',2)

% hold on
plot([0,0],[-1.1, 1.1],'k','LineWidth',3)
hold on
plot(disk_radius*sin(t)-sigma,disk_radius*cos(t),'r','LineWidth',2)
hold on
plot(real(eig_A),imag(eig_A),'ob')
hold off

axis equal

Expected_LLE_if_zero_mean = disk_radius-sigma
actual_LLE = max(real(eig(corrected_A)))

ylabel('Imag(\lambda_{i})','FontSize',18)
xlabel('Real(\lambda_{i})','FontSize',18)
% axis([-1.1 1.1 -1.1 1.1])

box off

% plot A (not shifted version)
figure(3)
set(gcf,'Position',[1983         325         560         420])
imagesc(A,[-2 2])
colorbar


% save_some_figs_to_folder_2('figs', 'girko_mu_0_sparse_0_sigma_0_b0p5',[1], [])

%% check if strongly connected
G = digraph(sparse(A));                  % build a directed graph from A
bins = conncomp(G,'Type','strong');      % strongly connected components
A_isStronglyConnected = (numel(unique(bins)) == 1)
nSCC = numnodes(condensation(G))

disk_radius = disk_radius  
