load('post_samples_real_new_resid_j24_pen_B_spline_cubic_2knots_aurora_n4000_sigma_j_sq_0.01_0.0002_nu4_long_init0.mat')
post_samples0 = post_samples;
load('post_samples_real_new_resid_j24_pen_B_spline_cubic_2knots_aurora_n4000_sigma_j_sq_0.01_0.0002_nu4_long_init1.mat')
post_samples1 = post_samples;

for i = 1:6
    subplot(3, 2, i)
    plot(post_samples0.eta(i,:))
    hold on
    plot(post_samples1.eta(i, :), 'r')
    title(['eta', num2str(i)])
end