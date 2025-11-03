%% if A is sparse, but mean in degree is held fixed, then spectral radius is constant

close all
clear
clc

set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 18);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

% dims to vary: n, density, meaan_in_out_degree, mu, 

n_vec = [250, 500]
n_plots = length(n_vec)

for i_n = 1:length(n_vec)

    n = n_vec(i_n);
    b = 2; % standard deviation
    mu = 0.0;
    
    mean_in_out_degree = 30; % desired mean number of connections in and out
    density = mean_in_out_degree/(n); % each neuron can make up to n-1 connections with other neurons

    disk_radius = sqrt(n)*b*sqrt(density)
    
    sigma = disk_radius;
    
    %% Generate A, Ginibre
    
    A = b*randn(n,n) + mu;
    
    A(rand(n,n)<(1-density)) = 0;
    
    % A(eye(size(A),'logical')) = 0; % no self connection
    
    %% Generate A, strongly connected
    
    
    %% check A stats
    
    actual_var_A = var(A(:))
    actual_stdev_A = std(A(:))
    actual_mean_A = mean(A(:))
    actual_mean_in_out_degree = sum(abs(A(:))>0)/n
    actual_density = sum(abs(A(:))>0)/n^2
    
    %%
    corrected_A = A-sigma*eye(size(A));
    corrected_var_A = var(A(:))
    corrected_stdev_A = std(A(:))
    corrected_mean_A = mean(A(:))
    
    eig_A = eig(corrected_A);
    
    figure(1)
    subplot(ceil(sqrt(n_plots)),ceil(sqrt(n_plots)),i_n)
    set(gcf,'Position',[1880         787         646         530])
    t = 0:.01:2*pi;
    % plot(sin(t)-sigma,cos(t),'g','LineWidth',2)
    
    % hold on
    plot([0,0],[-8, 8],'k','LineWidth',3)
    hold on
    plot(disk_radius*sin(t)-sigma,disk_radius*cos(t),'r','LineWidth',2)
    hold on
    plot(real(eig_A),imag(eig_A),'.b')
    hold off
    
    axis equal
    
    Expected_LLE_if_zero_mean = disk_radius-sigma
    actual_LLE = max(real(eig(corrected_A)))
    
    ylabel('Imag(\lambda_{i})','FontSize',18)
    xlabel('Real(\lambda_{i})','FontSize',18)
    title(['n = ' num2str(n)],'FontWeight','normal')
    % axis([-1.1 1.1 -1.1 1.1])
    
    box off

end

%%
save_some_figs_to_folder_2('RMT_review_figs', 'RMT_multipannel',[1], [])

% G = digraph(sparse(A));                  % build a directed graph from A
% bins = conncomp(G,'Type','strong');      % strongly connected components
% A_isStronglyConnected = (numel(unique(bins)) == 1)
% nSCC = numnodes(condensation(G))