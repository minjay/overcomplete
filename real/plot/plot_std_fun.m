% load data
load('data_EOF_regr_new.mat')

% non-stationary variance function
[~, ~, ~, X] = SCHA_regr(zeros(size(phi)), theta, phi, 3, 3, 4);
b_mat = X;

std_vec = exp(b_mat*post_samples.eta);

nu = 4;

for i = 1:100:size(post_samples.sigma_j_sq, 2)
    sigma_j  = sqrt(post_samples.sigma_j_sq(:, i));
    variance = plot_corr_fun_flex(2, sigma_j, 2, 3, nu, 1000, 'r');
    % rescale
    std_vec(:, i) = std_vec(:, i)*sqrt(variance)*1000;
end

figure
plot_pot(reshape(std_vec(:, i), size(phi)), phi, theta, 1000, max(abs(std_vec(:, i))))



