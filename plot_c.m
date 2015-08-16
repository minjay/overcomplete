function plot_c(c, post_samples_c)

c_mean = mean(post_samples_c, 2);
subplot(2, 1, 1)
y_min = min(min(c), min(c_mean));
y_max = max(max(c), max(c_mean));
plot(c)
axis([1 length(c) y_min y_max])
title('True values of c')
subplot(2, 1, 2)
plot(c_mean)
axis tight
axis([1 length(c) y_min y_max])
title('Estimated values of c')

end