% precompute the design matrix

load('theta_phi_R.mat')

B = 2;
j_min = 2;
j_max = 4;

theta_vec = theta(:);
phi_vec = phi(:);

% design matrix A
[Npix, ~, A_part_theta] = get_A_part_ss('theta', B, j_min, j_max, theta_vec*4, phi_vec);
[Npix, ~, A_part_phi] = get_A_part_ss('phi', B, j_min, j_max, theta_vec*4, phi_vec);


save('mat_A_part.mat', 'A_part_theta', 'A_part_phi', '-v7.3')


