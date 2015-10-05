load('theta_phi.mat')
load('post_samples_real.mat')

theta_vec = theta(:);
phi_vec = phi(:);

N = length(theta_vec);

eta_est = mean(post_samples.eta, 2);

r = 4;
mu = pi/(r+1)*(1:r);
lambda = pi/(r+1)*2.5/2;
b_mat = zeros(N, r+1);
b_mat(:, 1) = 1;
for i = 2:r+1
    b_mat(:, i) = exp(-(theta_vec*4-mu(i-1)).^2/2/lambda^2);
end

std_vec_est = exp(b_mat*eta_est);

n_samples = size(post_samples.eta, 2);
std_vec_post = zeros(N, n_samples);
for i = 1:n_samples
    std_vec_post(:, i) = exp(b_mat*post_samples.eta(:, i));
end

CI = quantile(std_vec_post', [0.025 0.975]);

[theta_sort, index] = sort(theta_vec);
lat_sort = theta_sort/pi*180;
plot(lat_sort, std_vec_est(index), 'b-', 'LineWidth', 1.5)
hold on
plot(lat_sort, CI(1, index), 'r--', 'LineWidth', 1.5)
plot(lat_sort, CI(2, index), 'r--', 'LineWidth', 1.5)
legend('Est. std function', '95% CI endpoints')
axis tight
