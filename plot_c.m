function plot_c(c, post_samples_c)

c_mean = mean(post_samples_c, 2);
Npix = [156 564 2148];
subplot(2, 2, 1)
plot(c(1:Npix(1)), c_mean(1:Npix(1)), '.')
axis square
xlabel('True')
ylabel('Fitted')
refline(1, 0);
subplot(2, 2, 2)
plot(c(Npix(1)+1:Npix(1)+Npix(2)), c_mean(Npix(1)+1:Npix(1)+Npix(2)), '.')
axis square
xlabel('True')
ylabel('Fitted')
refline(1, 0);
subplot(2, 2, 3)
plot(c(Npix(1)+Npix(2)+1:end), c_mean(Npix(1)+Npix(2)+1:end), '.')
axis square
xlabel('True')
ylabel('Fitted')
refline(1, 0);

end