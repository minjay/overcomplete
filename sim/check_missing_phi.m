clear

load('post_samples_missing_phi.mat')
load('ss.mat')

B = 2;
j_min = 2;
j_max = 3;

theta = 0:0.01:pi;

% non-stationary variance function
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta);

b_mat(:, 1) = 1;

eta_est = mean(post_samples.eta(:, 501:end), 2);
plot(theta, exp(b_mat*eta), 'LineWidth', 2)
hold on
plot(theta, exp(b_mat*eta_est), 'LineWidth', 2)
legend('True', 'Fitted')
set(gca, 'FontSize', 12)

sigma_j(2)
mean(sqrt(post_samples.sigma_j_sq(2, 501:end)))

c_est = mean(post_samples.c(:, 501:end), 2);

figure
st = 1;
for j = j_min:j_max 
    index_j = j-j_min+1;
    t = 2*floor(B^(j+1))+1;
    grid_points{index_j} = ss{degree_t==t};
    Npix = size(grid_points{index_j}, 1);
    x_grid = grid_points{index_j}(:, 1);
    y_grid = grid_points{index_j}(:, 2);
    z_grid = grid_points{index_j}(:, 3);
    [phi_grid, ~, ~] = cart2sph(x_grid, y_grid, z_grid);
    phi_grid(phi_grid<0) = phi_grid(phi_grid<0)+2*pi;
    c_j = c(st:st+Npix-1);
    c_est_j = c_est(st:st+Npix-1);
    index = find(phi_grid>3/2*pi);
    index2 = find(phi_grid<=3/2*pi);
    subplot(1, 2, index_j)
    plot(c_est_j(index2), c_j(index2), 'o')
    hold on
    plot(c_est_j(index), c_j(index), 'x')
    axis equal
    axis square
    legend('Data', 'No data', 'location', 'best')
    xlabel('Fitted c')
    ylabel('True c')
    title(['j = ', num2str(j)])
    st = st+Npix;
end


load('post_samples_missing_phi.mat')

theta = theta/4;
phi_rot = phi+pi/2;
[x_samples, y_samples] = pol2cart(phi_rot, theta/pi*180);

h = mypolar_grid([0 2*pi], [0 45]);
delete(h)
hold on
scatter(x_samples, y_samples, 20, 'bo')
