function pred_sim_needlet_kriging(seed, name)

load(['data_sim_', name, '.mat'])
load(['post_samples_', name, '_', num2str(seed), '.mat']) 

B = 2;
j_min = 2;
j_max = 4;

% sampling
N = length(Y);
n = length(index);
Y_samples = Y(index);

[x, y, z] = trans_coord(theta_vec, phi_vec);

% get inner prod matrix
r = zeros(N);
for j = 1:N
    for i = 1:j-1
        tmp = sum([x(i) y(i) z(i)].*[x(j) y(j) z(j)]);
        tmp = min(tmp, 1);
        tmp = max(tmp, -1);
        r(i, j) = tmp;
        r(j, i) = tmp;
    end
    r(j, j) = 1;
end

[r_vec, ~, ic] = unique(r(:));
n_r = length(r_vec);

% get cov vec
l_max = floor(B^(j_max+1));
sigma_j_est = mean(sqrt(post_samples.sigma_j_sq), 2);

Pl_mat = p_polynomial_value(n_r, l_max, r_vec);
cov_vec = zeros(n_r, 1);
for j = j_min:j_max
    l_min = ceil(B^(j-1));
    l_max = floor(B^(j+1));
    tmp = zeros(n_r, 1);
    for l = l_min:l_max
        tmp = tmp+(fun_b(l/B^j, B))^2*(2*l+1)/(4*pi)*Pl_mat(:, l+1);
    end
    cov_vec = cov_vec+nu*(sigma_j_est(j-j_min+1))^2/(nu-2)*tmp;
end

tau = mean(1./sqrt(post_samples.tau_sq_inv));
cov_mat = reshape(cov_vec(ic), N, N)+tau^2*eye(N);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(Y_samples, n, 1);

index_pred = setdiff(1:N, index);
SigmaP0 = cov_mat(:, index);
Y_pred_needlet_kriging = SigmaP0*tmp;

err = Y(index_pred)-Y_pred_needlet_kriging(index_pred);
mean(err.^2)

end
