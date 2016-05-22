load('theta_phi_R.mat')
load('mat_A.mat')
load('mat_A_part.mat')

A = A(361:16920, :);

plot_pot(reshape(A(:, 2), size(phi)), phi, theta, 1000, max(abs(A(:, 2))))

figure
plot_pot(reshape(A_part_theta(:, 2), size(phi)), phi, theta, 1000, max(abs(A_part_theta(:, 2))))
figure
plot_pot(reshape(A_part_phi(:, 2), size(phi)), phi, theta, 1000, max(abs(A_part_phi(:, 2))))
