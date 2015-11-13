load('data_sim.mat')

rng(1)

% sampling
N = length(Y);
n = 1e3;
index = randsample(N, n);
theta_samples = theta(index);
phi_samples = phi(index);
pot_samples = Y(index);

% non-stationary variance funcion
m = 4;
lambda_inv = 2/2.5;
b_mat = get_nonsta_var(m, lambda_inv, theta_samples);

beta_hat = [4.698266 1.532569 0.769803 0.832334 -0.860071 4.207303 17.918832 8.156696]; 
%beta_hat = [3.923432 2.04541 0.396874 1.506968 -0.581043 4.839481 19.966406 9.908376]; 
%beta_hat = [0.8937, 6.2331, -3.5112, 5.1321, 0.0806, 4.9487, 19.9509, 11.2629];
std_vec_est = exp(b_mat*reshape(beta_hat(1:m+1), m+1, 1));
plot(theta_samples, std_vec_est, '.')

n_r = 1e3;
r_vec = linspace(0, 2, n_r);
value_vec = zeros(n_r, 1);
nu = beta_hat(6);
a = beta_hat(7);
for i = 1:n_r
    value_vec(i) = Matern(r_vec(i), nu, a);
end
% convert chordal distance to great-circle distance
r_vec = asin(r_vec/2)*2;
% plot correlation function
plot(r_vec, value_vec, 'LineWidth', 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x, y, z] = trans_coord(theta_vec, phi_vec);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance funcion
b_mat = get_nonsta_var(m, lambda_inv, theta_vec);

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
cov_mat = get_cov_nonsta_Matern(beta, r, b_mat)+tau^2*eye(N);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(pot_samples, n, 1);

SigmaP0 = cov_mat(:, index);
Y_pred = SigmaP0*tmp;

save('Y_pred.mat', 'Y_pred')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L, T] = meshgrid(phi-pi, pi/2-theta);
[HX, HY] = sph2hammer(L, T);

% plot
figure
cmax = max(abs(Y_pred));
pcolor(HX, HY, reshape(Y_pred, size(phi_mat)));
shading interp
axis equal
axis tight
axis off
caxis([-cmax cmax])
