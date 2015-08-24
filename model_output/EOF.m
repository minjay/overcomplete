load('WHI_quad.mat')

all_Pot_N = all_Pot_N(1:720);
n = length(all_Pot_N);

% design matrix
X = zeros(n, numel(theta));
for i = 1:n
    X(i, :) = all_Pot_N{i}(:);
end

clear all_Pot_N

mean_map = mean(X, 1);
plot_pot(reshape(mean_map, size(phi)), phi, theta, 1000)

for i = 1:n
    X(i, :) = X(i, :)-mean_map;
end

[U, S, V] = svd(X, 'econ');

% find the 95% threshold
lambda = diag(S).^2;
var_exp = cumsum(lambda)./sum(lambda);
thres = find(var_exp>=0.95, 1, 'first');

% compute the residual fields
V_L = V(:, [1:thres]);
T_L = X*V_L;
X_L = T_L*V_L';
r = X-X_L;

cov_mat = corr(r);
theta_vec = theta(:);
phi_vec = phi(:);
[x, y, z] = sph2cart(phi_vec, pi/2-theta_vec, 1);
n = length(x);
dist_mat = zeros(n);
for i = 1:n
    if mod(i, 100)==0
        i
    end
    for j = i:n
        inner_product = sum([x(i), y(i), z(i)].*[x(j), y(j), z(j)]);
        inner_product = min(inner_product, 1);
        inner_product = max(inner_product, -1);
        dist_mat(i, j) = acos(inner_product);
        dist_mat(j, i) = acos(inner_product);
    end
end
        
dist_mat_sub = dist_mat(1:10:end, 1:10:end);
cov_mat_sub = cov_mat(1:10:end, 1:10:end);

k = 1;
figure
plot_pot(reshape(X(k, :), size(theta)), phi, theta, 1000)
figure
plot_pot(reshape(r(k, :), size(theta)), phi, theta, 1000)
figure
plot_pot(reshape(X(k, :)-r(k, :), size(theta)), phi, theta, 1000)

std_map = std(r, 0, 1);
plot_pot(reshape(std_map, size(theta)), phi, theta, 1000)

figure
r_field = reshape(r(k, :)./std_map, size(theta));
plot_pot(r_field, phi, theta, 1000)