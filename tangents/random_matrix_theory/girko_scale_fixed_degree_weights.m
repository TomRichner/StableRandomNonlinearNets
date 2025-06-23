close all
clear
clc

set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 18);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

n = 1000;
b = 1;
mu = 0.0;

mean_in_out_degree = 10; % desired mean number of connections in and out
density = mean_in_out_degree/(n-1); % each neuron can make up to n-1 connections with other neurons
sparse_fraction = 1-density;

disk_radius1 = b*sqrt(1-sparse_fraction);

% sigma = disk_radius1;
sigma = 0


% A = b*randn(n,n)./sqrt(n) + mu;
A = b*randn(n,n) + mu;
A(eye(size(A),'logical')) = 0; % no self connection

A(rand(n,n)<sparse_fraction) = 0;
A = A-sigma*eye(size(A));


var_A = var(A(:))
mean_A = mean(A(:))
eig_A = eig(A);



figure(1)
t = 0:.01:2*pi;
% plot(sin(t)-sigma,cos(t),'g','LineWidth',2)

% hold on
plot([0,0],[-1.1, 1.1],'k','LineWidth',3)
hold on
plot(disk_radius1*sin(t)-sigma,disk_radius1*cos(t),'r','LineWidth',2)
hold on
plot(real(eig_A),imag(eig_A),'ob')
hold off

axis equal

Expected_LLE = mu*n*sparse_fraction-sigma
ylabel('Imag(\lambda_{i})','FontSize',18)
xlabel('Real(\lambda_{i})','FontSize',18)
% axis([-1.1 1.1 -1.1 1.1])

box off

save_some_figs_to_folder_2('figs', 'girko_mu_0_sparse_0_sigma_0_b0p5',[1], [])