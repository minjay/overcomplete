load('data_EOF_regr_new.mat')

% plot field
figure
plot_pot(reshape(resid_all(1, :), size(phi)), phi, theta, 1000, max(abs(resid_all(1, :))))

% plot empirical std function
figure
emp_std_vec = std(resid_all, 0, 1);
theta_vec = theta(:);
plot(theta_vec, emp_std_vec, '.')