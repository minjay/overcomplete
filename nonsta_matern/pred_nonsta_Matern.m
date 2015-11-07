load('data.mat')

theta_vec = theta(:);
phi_vec = phi(:);
[x, y, z] = trans_coord(theta_vec, phi_vec);

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
    b_mat(:, i) = exp(-(theta_vec*4-mu(i-1)).^2/2/lambda^2);
end

beta = beta_hat(1:end-1);
tau = beta_hat(end);
cov_mat = get_cov(beta, r, b_mat)+eye(n)*tau^2;

Sigma00 = cov_mat(index, index);
tmp = Sigma00\pot_samples';

SigmaP0 = cov_mat(:, index);
Y_pred = SigmaP0*tmp;