function sim_data(nu, name)

rng(1)

alpha = 3;

tau = 0.1;

% the grid
res = 100;

theta = linspace(0, pi, res/2);
d_phi = 2*pi/res;
phi = linspace(0, 2*pi-d_phi, res);

[phi_mat, theta_mat] = meshgrid(phi, theta);
theta_vec = theta_mat(:);
phi_vec = phi_mat(:);

B = 2;
j_min = 2;
j_max = 3;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_vec, phi_vec);
[N, M] = size(A);

sigma_j = B.^(-alpha/2*(j_min:j_max));
sigma_j = sigma_j/sigma_j(1);

std_vec = exp(-(theta_vec-pi/2).^2/2/(pi/4)^2);
DA = zeros(N, M);
for i = 1:N
    DA(i, :) = std_vec(i)*A(i, :);
end

fj_sq = zeros(M, 1);
c = zeros(M, 1);
st = 1;
for j = j_min:j_max
    index_j = j-j_min+1;
    range = st:st+Npix(index_j)-1;
    fj_sq(range) = sigma_j(index_j)^2*ones(Npix(index_j), 1);
    c(range) = sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
    st = st+Npix(index_j);
end

Y = DA*c+tau*randn(N, 1);

save(['data_sim_', name, '.mat'], 'theta', 'phi', 'theta_mat', 'phi_mat', 'theta_vec',...
    'phi_vec', 'Y', 'nu', 'alpha', 'tau', 'sigma_j', 'fj_sq')

end