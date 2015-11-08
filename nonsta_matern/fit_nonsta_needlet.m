parpool

load('data.mat')

rng(1)

% sampling
N = 1e3;
[pot_samples, theta_samples, phi_samples, index] = sampling_data(resid,...
    theta, phi, N, 0);

[x, y, z] = trans_coord(theta_samples*4, phi_samples);

n = length(x);
r = zeros(n);
for j = 1:n
    for i = 1:j
        s = [x(i); y(i); z(i)];
        t = [x(j); y(j); z(j)];
        inner_prod = dot(s, t);
        inner_prod = max(inner_prod, -1);
        inner_prod = min(inner_prod, 1);
        r(i, j) = inner_prod;
    end
end

% non-stationary variance funcion
m = 4;
mu = pi/(m+1)*(1:m);
lambda = pi/(m+1)*2.5/2;
b_mat = zeros(n, m+1);
b_mat(:, 1) = 1;
for i = 2:m+1
    b_mat(:, i) = exp(-(theta_samples*4-mu(i-1)).^2/2/lambda^2);
end

B = 2;
jj = 2;
beta_all = [zeros(5, 1); 0.01];
f_value = negloglik_needlet(beta_all, r, b_mat, pot_samples', B, jj);

if ~paral
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point',...
        'Hessian', 'lbfgs', 'Display', 'iter', 'TolX', 1e-6);
else
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point',...
        'Hessian', 'lbfgs', 'Display', 'iter', 'TolX', 1e-6,...
        'UseParallel', 'always');
end

beta_init = [zeros(5, 1); 2; 1; 0.01];
negloglik1 = @(beta_all) negloglik_needlet(beta_all, r, b_mat, pot_samples');

lb = [-Inf -Inf -Inf -Inf -Inf 0 0 0];
ub = [Inf Inf Inf Inf Inf 5 Inf Inf];

[beta_hat, f_min] = fmincon(negloglik1, beta_init, [], [], [], [], lb, ub, [], options);

beta_hat = [3.1530, 5.6248, -2.9758, 4.3551, -0.0394, 2.7891, 6.0510, 10.5671];
std_vec_est = exp(b_mat*beta_hat(1:5)');
plot(theta_samples, std_vec_est, '.')

r_vec = linspace(0, 2, 1000);
value = zeros(length(r_vec), 1);
nu = beta_hat(6);
a = beta_hat(7);
tau = beta_hat(8);
for i = 1:length(r_vec)
    value(i) = Matern(r_vec(i), nu, a);
end
r_vec = asin(r_vec/2)*2;
plot(r_vec, value, 'LineWidth', 2)