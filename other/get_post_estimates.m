load('post_samples.mat')

rng(1)

% the grid
B = 2;
% change Nside to 8 or 16
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
theta = zeros(N, 1);
phi = zeros(N, 1);
for i = 1:N
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end

% perturbation
theta = theta+randn(N, 1)*pi/10;
theta(theta<0) = theta(theta<0)+pi;
theta(theta>pi) = theta(theta>pi)-pi;
phi = phi+randn(N, 1)*2*pi/10;
phi(phi<0) = phi(phi<0)+2*pi;
phi(phi>2*pi) = phi(phi>2*pi)-2*pi;

j_min = 2;
j_max = 4;

% design matrix A
[Npix, ~, ~] = get_A_ss(B, j_min, j_max, theta, phi);

eta_est = mean(post_samples.eta, 2);

r = 4;
mu = pi/(r+1)*(1:r);
lambda = pi/(r+1)*2.5/2;
b_mat = zeros(N, r+1);
b_mat(:, 1) = 1;
for i = 2:r+1
    b_mat(:, i) = exp(-(theta-mu(i-1)).^2/2/lambda^2);
end

std_vec_true = exp(b_mat*eta);
std_vec_est = exp(b_mat*eta_est);

n_samples = size(post_samples.eta, 2);
std_vec_post = zeros(N, n_samples);
for i = 1:n_samples
    std_vec_post(:, i) = exp(b_mat*post_samples.eta(:, i));
end

CI = quantile(std_vec_post', [0.05 0.95]);

[theta_sort, index] = sort(theta);
plot(theta_sort, std_vec_true(index), 'LineWidth', 1.5)
hold on
plot(theta_sort, std_vec_est(index), 'r--', 'LineWidth', 1.5)
plot(theta_sort, CI(1, index), 'g-.', 'LineWidth', 1.5)
plot(theta_sort, CI(2, index), 'g-.', 'LineWidth', 1.5)
legend('True std function', 'Est. std function', '95% CI endpoints')
axis tight

% plot c
c_mean = mean(post_samples.c, 2);
loc_x = [0.1 0.4 0.7];
loc_y = [0.1 0.1 0.1];
st = 1;
for j = j_min:j_max
    index_j = j-j_min+1;
    range = st:st+Npix(index_j)-1;
    subplot('Position', [loc_x(index_j) loc_y(index_j) 0.2 0.8])
    plot(c(range), c_mean(range), '.')
    xlabel('True c')
    ylabel('Est. c')
    title(['j = ', num2str(j)])
    axis tight
    axis square
    hline = refline(1, 0);
    set(hline, 'color', 'red', 'LineStyle', '--', 'LineWidth', 1.5);
    st = st+Npix(index_j);
end

% boxplot
boxplot([post_samples.eta' 1./sqrt(post_samples.tau_sq_inv)'])
% add labels
set(gca, 'XTickLabel', {' '})
text(1, -2, '\eta_0')
text(2, -2, '\eta_1')
text(3, -2, '\eta_2')
text(4, -2, '\eta_3')
text(5, -2, '\eta_4')
text(6, -2, '\tau')