clear

% load fitted result
load('post_samples_real_exp3.mat')

% load data
load('data_EOF_regr_new.mat')
resid = resid_all(1, :);

figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.05], [0.05 0.05], [0.05 0.02]);

for i = 1:4
    subplot(2, 4, i)
    plot(post_samples.eta(i, 1:3000))
    hold on
    y_lim = get(gca, 'ylim');
    plot([2000 2000], y_lim, 'LineWidth', 2)
    axis tight
    title(['\eta_', num2str(i-1)])
end
subplot(2, 4, 5)
plot(post_samples.sigma_j_sq(2, 1:3000))
hold on
y_lim = get(gca, 'ylim');
plot([2000 2000], y_lim, 'LineWidth', 2)
axis tight
title('\sigma_3^2')
subplot(2, 4, 6)
plot(post_samples.sigma_j_sq(3, 1:3000))
hold on
y_lim = get(gca, 'ylim');
plot([2000 2000], y_lim, 'LineWidth', 2)
axis tight
title('\sigma_4^2')
subplot(2, 4, 7)
plot(post_samples.tau_sq_inv(1:3000))
hold on
y_lim = get(gca, 'ylim');
plot([2000 2000], y_lim, 'LineWidth', 2)
axis tight
title('1/\tau^2')
subplot(2, 4, 8)
plot((post_samples.acc_times_all(1:2:6000)+post_samples.acc_times_all(2:2:6000))/2)
hold on
y_lim = get(gca, 'ylim');
plot([2000 2000], y_lim, 'LineWidth', 2)
axis tight
title('Acceptance rate')

% scatter plot
[h, ax] = plotmatrix([post_samples.eta(:, 2001:3000)' post_samples.sigma_j_sq(2:3, 2001:3000)'...
    post_samples.tau_sq_inv(2001:3000)']);
labels = {'\eta_0', '\eta_1', '\eta_2', '\eta_3', '\sigma_3^2', '\sigma_4^2',...
    '1/\tau^2'};
% make labels
for i = 1:length(ax)
    xlabel(ax(end, i), labels{i})
    ylabel(ax(i, 1), labels{i})
end

print -painters -depsc scatter_plot_real.eps

subplot = @(m,n,p) subtightplot (m, n, p, [0.125 0.075], [0.075 0.05], [0.075 0.02]);

% autocorr
for i = 1:4
    subplot(2, 4, i)
    autocorr(post_samples.eta(i, 2001:3000))
    axis tight
    title(['\eta_', num2str(i-1)])
end
subplot(2, 4, 5)
autocorr(post_samples.sigma_j_sq(2, 2001:3000))
axis tight
title('\sigma_3^2')
subplot(2, 4, 6)
autocorr(post_samples.sigma_j_sq(3, 2001:3000))
axis tight
title('\sigma_4^2')
subplot(2, 4, 7)
autocorr(post_samples.tau_sq_inv(2001:3000))
axis tight
title('1/\tau^2')
