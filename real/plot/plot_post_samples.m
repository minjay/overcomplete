% load post_samples

% plot post_samples
figure
for i = 1:16
    subplot(4, 4, i)
    plot(post_samples.eta(i, :))
    title(['eta ', num2str(i)])
end

figure
plot(post_samples.sigma_j_sq(2, :))