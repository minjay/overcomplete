function sim_data(nu, name)
% function to simulate data on a perturbed HEALPix grid

rng(1)

alpha = 3;

tau = 0.1;

% the grid
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
theta0 = zeros(N, 1);
phi0 = zeros(N, 1);
for i = 1:N
    theta0(i) = tp{i}(1);
    phi0(i) = tp{i}(2);
end

% perturbation
theta = theta0+randn(N, 1)*pi/10;
theta(theta<0) = theta(theta<0)+pi;
theta(theta>pi) = theta(theta>pi)-pi;
phi = phi0+randn(N, 1)*2*pi/10;
phi(phi<0) = phi(phi<0)+2*pi;
phi(phi>2*pi) = phi(phi>2*pi)-2*pi;

B = 2;
j_min = 2;
j_max = 3;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta, phi);
[N, M] = size(A);

sigma_j = B.^(-alpha/2*(j_min:j_max));
sigma_j = sigma_j/sigma_j(1);

% non-stationary variance function
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta);

b_mat(:, 1) = 1;

r = size(b_mat, 2)-1;

rng(2)
eta = randn(r+1, 1);
std_vec = exp(b_mat*eta);
DA = zeros(N, M);
for i = 1:N
    DA(i, :) = std_vec(i)*A(i, :);
end

fj_sq = zeros(M, 1);
c = zeros(M, 1);
st = 1;
for j = j_min:j_max
    index_j = j-j_min+1;
    range = st:st+Npix(index_j)-1;
    fj_sq(range) = sigma_j(index_j)^2*ones(Npix(index_j), 1);
    c(range) = sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
    st = st+Npix(index_j);
end

Y = DA*c+tau*randn(N, 1);

save(['data_sim_', name, '.mat'], 'theta', 'phi', 'theta',...
    'phi', 'Y', 'nu', 'alpha', 'tau', 'sigma_j', 'fj_sq')

end