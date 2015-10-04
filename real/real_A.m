B = 2;
j_min = 2;
j_max = 4;
nu = 4;
alpha = 4;

res = 500;

phi = linspace(0, 2*pi, res);
theta = linspace(0, pi, res/2);
[theta_mat, phi_mat] = meshgrid(theta, phi);
theta_vec = theta_mat(:);
phi_vec = phi_mat(:);

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_vec, phi_vec);

save('mat_A.mat', 'A', '-v7.3')


