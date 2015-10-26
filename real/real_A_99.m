% precompute the design matrix

load('theta_phi.mat')

B = 2.5;
j_min = 2;
j_max = 3;
nu = 4;
alpha = 4;

theta_vec = theta(:);
phi_vec = phi(:);

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_vec*4, phi_vec);

save('mat_A_99.mat', 'A', '-v7.3')


