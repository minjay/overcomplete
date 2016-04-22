% load data
load('data_EOF_regr_new.mat')
% load pre-computed design matrix
load('mat_A.mat')

% set seed
rng(1)

% parameter setting
% d.f.
nu = 4;
% j's
j_min = 2;
j_max = 4;

% convert theta to vector
theta_vec = theta(:);

% simulate
% non-stationary variance function
knots = [0 0 0 0 40/180 80/180 120/180 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_vec*4);

% discard burn-in period
burn_in = 500;
post_samples_eta = post_samples.eta(:, burn_in+1:end);
post_samples_sigma_j = sqrt(post_samples.sigma_j_sq(:, burn_in+1:end));
post_samples_tau = 1./sqrt(post_samples.tau_sq_inv(burn_in+1:end));

% num of posterior samples
n_sample = size(post_samples_eta, 2);

% posterior samples of std_vec
std_vec = exp(b_mat*post_samples_eta);

% subset A
[N, M] = size(A);
A = A(361:(N-360), :);
[N, M] = size(A);

Y = zeros(N, n_sample);

for i = 1:n_sample
    i
    DA = zeros(N, M);
    for j = 1:N
        DA(j, :) = std_vec(j, i)*A(j, :);
    end
    c = zeros(M, 1);
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range) = post_samples_sigma_j(index_j, i)*trnd(nu, Npix(index_j), 1);
        st = st+Npix(index_j);
    end
    Y(:, i) = DA*c+post_samples_tau(i)*randn(N, 1);
end

% rescale
Y = Y*1e3;

% save
save('uncond_sim.mat', 'Y')