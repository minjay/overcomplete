load('data.mat')

rng('shuffle')
N = 500;
[pot_samples, theta_samples, phi_samples, index] = sampling_data(resid,...
    theta, phi, N, 1);

% specify parameters
B = 2;
j_min = 2;
j_max = 4;
nu = 4;
alpha = 4;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_samples*4, phi_samples);
M = size(A, 2);

% non-stationary variance funcion
r = 4;
mu = pi/(r+1)*(1:r);
lambda = pi/(r+1)*2.5/2;
b_mat = zeros(N, r+1);
b_mat(:, 1) = 1;
for i = 2:r+1
    b_mat(:, i) = exp(-(theta_samples*4-mu(i-1)).^2/2/lambda^2);
end

T = 1000;
Y = zeros(N, T);
for t = 1:T
    if mod(t, 10)==0
        t
    end
    c = post_samples.c(:, t);
    eta = post_samples.eta(:, t);
    tau_sq = 1/post_samples.tau_sq_inv(t);
    std_vec = exp(b_mat*eta);
    DA = zeros(N, M);
    for i = 1:N
        DA(i, :) = std_vec(i)*A(i, :);
    end
    DAc = DA*c;
    Y(:, t) = mvnrnd(DAc, tau_sq*eye(N))';
end

Y_pred = mean(Y, 2)*1e3;
plot(pot_samples, Y_pred, '.')
refline(1, 0)
axis square