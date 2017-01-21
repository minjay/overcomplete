% script to summarize the estimated parameters

clear
clc

load('post_samples_real_exp3.mat')

range = 2001:3000;
eta = post_samples.eta(:, range);
sigma_j = sqrt(post_samples.sigma_j_sq(:, range));
tau = 1./sqrt(post_samples.tau_sq_inv(range));

disp('eta')
mean(eta, 2)
% use default normalization
std(eta, 0, 2)

disp('sigma_j')
mean(sigma_j, 2)
std(sigma_j, 0, 2)

disp('tau')
mean(tau, 2)
std(tau, 0, 2)