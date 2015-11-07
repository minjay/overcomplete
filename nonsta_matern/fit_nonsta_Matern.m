parpool(8)

load('data.mat')

rng(1)

% sampling
N = 1e3;
[pot_samples, theta_samples, phi_samples, index] = sampling_data(resid,...
    theta, phi, N, 0);

[x, y, z] = trans_coord(theta_samples, phi_samples);

n = length(x);
r = zeros(n);
for j = 1:n
    for i = 1:j
        s = [x(i); y(i); z(i)];
        t = [x(j); y(j); z(j)];
        r(i ,j) = norm(s-t);
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

beta_all = [zeros(5, 1); 2; 1; 0.01];
f_value = negloglik(beta_all, r, b_mat, pot_samples');

if ~paral
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point',...
        'Hessian', 'lbfgs', 'Display', 'iter', 'TolX', 1e-6);
else
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point',...
        'Hessian', 'lbfgs', 'Display', 'iter', 'TolX', 1e-6,...
        'UseParallel', 'always');
end

beta_init = [zeros(5, 1); 2; 1; 0.01];
negloglik1 = @(beta_all) negloglik(beta_all, r, b_mat, pot_samples');

lb = [-Inf -Inf -Inf -Inf -Inf 0 0 0];
ub = [Inf Inf Inf Inf Inf 5 Inf Inf];

[beta_hat, f_min] = fmincon(negloglik1, beta_init, [], [], [], [], lb, ub, [], options);

std_vec_est = exp(b_mat*beta_hat(1:5));
plot(theta_samples, std_vec_est, '.')