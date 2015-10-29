clear

res = 500;

B = 2;
j_min = 2;
j_max = 4;
len_j = j_max-j_min+1;

% grid points
theta = linspace(0, pi, res/2);
phi = linspace(0, 2*pi, res);

[phi_mat, theta_mat] = meshgrid(phi, theta);
theta_vec = theta_mat(:);
phi_vec = phi_mat(:);

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_vec, phi_vec);

save('mat_A_sim.mat', 'Npix', 'A', '-v7.3')