addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

clear
rng(1)

nu = 2.1;
alpha = 3;

tau = 0.1;

% the grid
B = 2;
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
theta = zeros(N, 1);
phi = zeros(N, 1);
for i = 1:N
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end

% perturbation
theta = theta+randn(N, 1)*pi/10;
theta(theta<0) = theta(theta<0)+pi;
theta(theta>pi) = theta(theta>pi)-pi;
phi = phi+randn(N, 1)*2*pi/10;
phi(phi<0) = phi(phi<0)+2*pi;
phi(phi>2*pi) = phi(phi>2*pi)-2*pi;

j_min = 2;
j_max = 3;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta, phi);
M = size(A, 2); 

sigma_j = B.^(-alpha/2*(j_min:j_max));
sigma_j = sigma_j/sigma_j(1);

% non-stationary variance function
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta);

b_mat(:, 1) = 1;

r = size(b_mat, 2)-1;

rng(3)
eta = randn(r+1, 1);
eta(1) = 2.0;
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

beta_init = [zeros(1, r+1) 0.1^2 1e-2];
negloglik1 = @(beta_all) negloglik_Gaussian_needlet(beta_all, b_mat, Y, Npix, A);

lb = [-10*ones(1, r+1) 0 1e-3];
ub = [10*ones(1, r+1) 1 Inf];

[beta_hat, f_min] = Gaussian_needlet_fit(negloglik1, beta_init, lb, ub, true);

save('beta_hat_init_value.mat', 'beta_hat')

delete(gcp)