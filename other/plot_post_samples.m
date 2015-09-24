load('post_samples.mat')

r = size(post_samples.eta, 1)-1;

for i = 1:r+1
    figure
    subplot('Position', [0.1 0.1 0.35 0.8])
    plot(post_samples.eta(i, :))
    subplot('Position', [0.6 0.1 0.35 0.8])
    autocorr(post_samples.eta(i, :))
    print(['eta', num2str(i-1)], '-dpng')
end

figure
subplot('Position', [0.1 0.1 0.35 0.8])
plot(post_samples.tau_sq_inv)
subplot('Position', [0.6 0.1 0.35 0.8])
autocorr(post_samples.tau_sq_inv)
print('tau', '-dpng')
