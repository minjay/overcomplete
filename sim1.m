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

j_min = 2;
j_max = 4;
j_len = j_max-j_min+1;
n_dist = 1e3;

[Npix, A] = get_A_ss(B, j_min, j_max, theta, phi, n_dist);
[N, M] = size(A);
ATA = A'*A;

sigma_j = sigma*B.^(-alpha*(j_min:j_max));
fj_sq = zeros(M, 1);
c = zeros(M, 1);
for j = j_min:j_max
    index = j-j_min+1;
    range = sum(Npix(1:index))-Npix(index)+1:sum(Npix(1:index));
    fj_sq(range) = B^(-2*alpha*j)*ones(Npix(index), 1);
    c(range) = sigma_j(index)*trnd(nu, Npix(index), 1);
end

Y = A*c+randn(N, 1)*tau;

ATY = A'*Y;

sigma0 = 1;
sigma0_sq = sigma0^2;
tau0 = 0.01;
tau0_sq_inv = 1/tau0^2;
V0_inv = ones(M, 1); 
T = 30000;
burn_in = 25000;
thin = 5;
post_samples = Gibbs_sampler(A, ATA, Y, ATY, fj_sq, nu, sigma0_sq, tau0_sq_inv, V0_inv, T, burn_in, thin);
    