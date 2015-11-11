clear
rng(1)

nu = 4;
alpha = 4;

tau = 0.1;

% the grid
B = 2;
Nside = 16;
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
j_max = 4;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta, phi);
M = size(A, 2); 

sigma_j = B.^(-alpha/2*(j_min:j_max));

% non-stationary variance funcion
r = 4;
mu = pi/(r+1)*(1:r);
lambda = pi/(r+1)*2.5/2;
b_mat = zeros(N, r+1);
b_mat(:, 1) = 1;
for i = 2:r+1
    b_mat(:, i) = exp(-(theta-mu(i-1)).^2/2/lambda^2);
end

rng(2)
eta = [1.5; randn(r, 1)];
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
    fj_sq(range) = B^(-alpha*j)*ones(Npix(index_j), 1);
    c(range) = sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
    st = st+Npix(index_j);
end

Y = DA*c+randn(N, 1)*tau;

save('data_sim.mat', 'theta', 'phi', 'Y')