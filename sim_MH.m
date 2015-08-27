clear

nu = 3;

alpha = 4;
tau = 0.1;

% the grid
B = 2;
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
n = length(tp);
theta = zeros(n, 1);
phi = zeros(n, 1);
for i = 1:n
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end

rng(1)
theta = theta+randn(n, 1)*pi/10;
theta(theta<0) = theta(theta<0)+pi;
theta(theta>pi) = theta(theta>pi)-pi;
phi = phi+randn(n, 1)*2*pi/10;
phi(phi<0) = phi(phi<0)+2*pi;
phi(phi>2*pi) = phi(phi>2*pi)-2*pi;

j_min = 2;
j_max = 3;
j_len = j_max-j_min+1;
n_dist = 1e3;

[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta, phi, n_dist);
[N, M] = size(A);

sigma_j = B.^(-alpha/2*(j_min:j_max));
fj_sq = zeros(M, 1);
c = zeros(M, 1);
for j = j_min:j_max
    index = j-j_min+1;
    range = sum(Npix(1:index))-Npix(index)+1:sum(Npix(1:index));
    fj_sq(range) = B^(-alpha*j)*ones(Npix(index), 1);
    c(range) = sigma_j(index)*trnd(nu, Npix(index), 1);
end

r = 1;
mu = pi/(r+1)*(1:r);
lambda = pi/(r+1)*2.5/2;
b_mat = zeros(r+1, N);
b_mat(1, :) = 1;
for i = 2:r+1
    b_mat(i, :) = normpdf(theta, mu(i-1), lambda);
end

eta = randn(r+1, 1);

Y = diag(exp(b_mat'*eta))*A*c+randn(N, 1)*tau;

tau0 = 0.01;
tau0_sq_inv = 1/tau0^2;
V0_inv = ones(M, 1); 
eta0 = eta;
tau_eta_sq = 1;
sigma_eta_sq = 0.075;
T = 150000;
n_report = 100;
burn_in = 50000;
thin = 100;

tic
post_samples = Gibbs_sampler_MH(A, Y, b_mat, fj_sq, nu, tau0_sq_inv,...
    V0_inv, eta0, tau_eta_sq, sigma_eta_sq, T, burn_in, thin, n_report);
toc
