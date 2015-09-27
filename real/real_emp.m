load('svd.mat')
load('theta_phi.mat')

n_d = 7;
n_thin = 20;

corr_mat = corr(r(1:720*n_d, 1:n_thin:end));

phi_vec = phi(:);
theta_vec = theta(:);
phi_vec = phi_vec(1:n_thin:end);
theta_vec = theta_vec(1:n_thin:end);
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

[X_MED, Y_MED, Y_LOW, Y_HIGH] = binned_plot(dist_vec, corr_vec);
figure
plot(X_MED, Y_MED, '.')

hold on
plot_corr_fun(2, 2, 4, 4, 4, 10000, 'r-')
line([0 pi], [0 0], 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5)
axis tight