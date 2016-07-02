clear

load('post_samples_real_new_resid_j24_pen_B_spline_cubic_2knots_aurora_n4000_nu3_long_init0.mat')
post_samples1 = post_samples;
load('post_samples_real_new_resid_j24_pen_B_spline_cubic_2knots_aurora_n4000_nu4_long_init0.mat')
post_samples2 = post_samples;

for i = 1:6
    subplot(3, 3, i)
    plot(post_samples1.eta(i, :))
    hold on
    plot(post_samples2.eta(i, :))
end

for i = 1:2
    subplot(3, 3, 6+i)
    plot(post_samples1.sigma_j_sq(i+1, :))
    hold on
    plot(post_samples2.sigma_j_sq(i+1, :))
end

subplot(3, 3, 9)
plot(post_samples1.tau_sq_inv)
hold on
plot(post_samples2.tau_sq_inv)
legend('nu=3', 'nu=4', 'Location', 'best')

