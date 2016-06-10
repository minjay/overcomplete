clear

load('data_EOF_regr_new.mat')
r = resid_all;

n_thin = 20;

phi_vec = phi(:);
theta_vec = theta(:);

phi_vec = phi_vec(1:n_thin:end);
theta_vec = theta_vec(1:n_thin:end);
corr_mat = corr(r(:, 1:n_thin:end));

[x, y, z] = sph2cart(phi_vec, pi/2-theta_vec*4, 1);

n = length(x);
dist_vec = zeros((n+1)*n/2, 1);
corr_vec = zeros((n+1)*n/2, 1);
index = 0;

for i = 1:n
    for j = i:n
        index = index+1;
        inner_prod = sum([x(i) y(i) z(i)].*[x(j) y(j) z(j)]);
        inner_prod = min(inner_prod, 1);
        inner_prod = max(inner_prod, -1);
        dist_vec(index) = acos(inner_prod);
        corr_vec(index) = corr_mat(i, j);
    end
end

plot(dist_vec, corr_vec, '.')
hold on

%%% plot fitted corr function by nonGau-need
B = 2;
j_min = 2;
j_max = 4;

nu = 4;
sigma_j = sqrt([1 0.01 0.0002]);

plot_corr_fun_flex(2, sigma_j, 2, 4, nu, 1000, 'r');

%%% plot fitted corr function by Gau-need
load('beta_hat_needlet.mat')
sigma_j = sqrt([1 beta_hat(7:8)]);
plot_corr_fun_flex(2, sigma_j, 2, 4, nu, 1000, 'g');

%%% plot fitted corr function by Matern
load('beta_hat.mat')
m = 6;
n_r = 1e3;
r_vec = linspace(0, 2, n_r);
corr_vec = zeros(n_r, 1);
nu = beta_hat(m+1);
a = beta_hat(m+2);
for i = 1:n_r
    corr_vec(i) = Matern(r_vec(i), nu, a);
end
% convert chordal distance to great-circle distance
r_vec = asin(r_vec/2)*2;
% plot correlation function
plot(r_vec, corr_vec, 'c', 'LineWidth', 2)

lg = legend('Empirical corr function', 'Est. corr function (nonGau-need)',...
    'Est. corr function (Gau-need)', 'Est. corr function (Gau-Matern)');
set(lg, 'FontSize', 12)
set(gca, 'FontSize', 12)
xlabel('Great-circle distance (rad)')
ylabel('Correlation')
axis tight