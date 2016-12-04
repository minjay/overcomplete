function [eta_est, sigma_j_est, tau_est] = sim_AM_sigma_rep(R, seed, eta, sigma_j, tau)
% simulation study
% generate non-stationary std function using random seed eta_seed

nu = 4;

% the grid
B = 2;
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
theta0 = zeros(N, 1);
phi0 = zeros(N, 1);
for i = 1:N
    theta0(i) = tp{i}(1);
    phi0(i) = tp{i}(2);
end

j_min = 2;
j_max = 3;

r = 4;

eta_est = zeros(r+1, R);
sigma_j_est = zeros(j_max-j_min+1, R);
tau_est = zeros(1, R);

rng(seed)
for rep = 1:R
    rep
    % perturbation
    theta = theta0+randn(N, 1)*pi/10;
    theta(theta<0) = theta(theta<0)+pi;
    theta(theta>pi) = theta(theta>pi)-pi;
    phi = phi0+randn(N, 1)*2*pi/10;
    phi(phi<0) = phi(phi<0)+2*pi;
    phi(phi>2*pi) = phi(phi>2*pi)-2*pi;

    % design matrix A
    [Npix, ~, A] = get_A_ss(B, j_min, j_max, theta, phi);
    M = size(A, 2); 

    % non-stationary variance function
    [b_mat, ~] = bspline_basismatrix(4, knots, theta);

    b_mat(:, 1) = 1;

    std_vec = exp(b_mat*eta);
    DA = zeros(N, M);
    for i = 1:N
        DA(i, :) = std_vec(i)*A(i, :);
    end

    c = zeros(M, 1);
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range) = sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
        st = st+Npix(index_j);
    end

    Y = DA*c+randn(N, 1)*tau;
    
    % get init values for MCMC
    beta_init = [zeros(1, r+1) 0.1^2 1e-2];
    negloglik1 = @(beta_all) negloglik_Gaussian_needlet(beta_all, b_mat, Y, Npix, A);

    lb = [-10*ones(1, r+1) 0 1e-3];
    ub = [10*ones(1, r+1) 1 Inf];

    [beta_hat, f_min] = Gaussian_needlet_fit(negloglik1, beta_init, lb, ub, false);

    % init
    % c
    c_init = zeros(M, 1);
    % V
    V_inv_init = ones(M, 1); 
    % sigma_j_sq
    sigma_j_sq_init = beta_hat(end-1);
    % eta
    eta_init = beta_hat(1:r+1)';
    % pri_sig of eta_0
    tau_sigma_sq = 1e2;
    % pri_sig of eta
    tau_eta_sq = 1e2;
    % tau
    tau_init = beta_hat(end);
    tau_sq_inv_init = 1/tau_init^2;
    % tuning parameters
    mu_init = zeros(r+1, 1);
    Sigma_init = eye(r+1);
    lambda = 0.05;
    % the number of MCMC iterations
    T = 4e5;
    % the length of the burn-in period
    burn_in = 2e5;
    % the length of the thinning interval
    thin = 200;
    % the length of the interval to report progress
    n_report = 1e3;

    model = struct('A', A, 'b_mat', b_mat, 'nu', nu);

    data = struct('Y', Y, 'Npix', Npix);

    params = struct('c', c_init, 'V', V_inv_init, 'sigma_j_sq', sigma_j_sq_init,...
        'eta', eta_init, 'tau_sigma_sq', tau_sigma_sq, 'tau_eta_sq', tau_eta_sq,...
        'tau', tau_sq_inv_init);

    tuning = struct('mu', mu_init, 'Sigma', Sigma_init, 'lambda', lambda);

    options = struct('T', T, 'burn_in', burn_in, 'thin', thin, 'n_report', n_report);

    post_samples = Gibbs_sampler_AM_rep_inter(model, data, params, tuning, options);
    
    eta_est(:, rep) = mean(post_samples.eta, 2);
    sigma_j_est(:, rep) = mean(sqrt(post_samples.sigma_j_sq), 2);
    tau_est(rep) = mean(1./sqrt(post_samples.tau_sq_inv));
    
end

end