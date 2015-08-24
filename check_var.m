clear

nu = 2;

alpha = 3;
sigma = 10;
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

j_min = 0;
j_max = 4;
j_len = j_max-j_min+1;
n_dist = 1e3;

[Npix, grid_points, A] = get_A_ss(B, j_min, j_max, theta, phi, n_dist);
[N, M] = size(A);

A_sq = A.*A;
sigma_j = sigma*B.^(-alpha*(j_min:j_max));
var_vec = zeros(N, 1);
for j = 4
    index = j-j_min+1;
    range = sum(Npix(1:index))-Npix(index)+1:sum(Npix(1:index));
    var_vec = var_vec+sum(A_sq(:, range), 2);
end

j = 2;
psi = zeros(N, Npix(j+1));
for k = 1:Npix(j+1)
    [~, psi_tmp] = get_psi_ss(B, j, k, theta, phi);
    psi(:, k) = psi_tmp;
end
var_vec2 = sum(psi.^2, 2);