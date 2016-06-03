% R is the radius of the ionosphere

load('theta_phi_R.mat')
load('deriv_B_spline.mat')
load('mat_A.mat')
load('mat_A_part.mat')
load('post_samples_real_new_resid_j24_pen_B_spline_cubic_2knots_aurora_n4000_sigma_j_sq_0.01_0.0002_nu4_long_init0.mat')
load('beta_hat_needlet.mat')

% set seed
rng(1)

% parameter setting
% j's
j_min = 2;
j_max = 4;

% convert theta to vector
theta_vec = theta(:);

% simulate
% non-stationary variance function
knots = [0 0 0 0 40/180 80/180 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_vec*4);
b_mat(:, 1) = 1;

m = size(b_mat, 2);

% posterior samples of std_vec
eta = beta_hat(1:m)';
std_vec = exp(b_mat*eta);

sigma_j = [1 sqrt(beta_hat(m+1:end-1))];

% subset A
[N, M] = size(A);
A = A(361:(N-360), :);
[N, M] = size(A);

b_mat_theta = kron(bS, ones(360, 1));
std_vec_theta = exp(b_mat*eta).*(b_mat_theta*eta);

T = 1e3;
neg_Y_theta_Gau_need = zeros(N, T);
neg_Y_phi_Gau_need = zeros(N, T);
Y_Gau_need = zeros(N, T);
DA_theta = zeros(N, M);
DA_phi = zeros(N, M);
D_thetaA = zeros(N, M);
DA = zeros(N, M);

for t = 1:T
    
    for j = 1:N
        DA_theta(j, :) = std_vec(j)*A_part_theta(j, :);
        DA_phi(j, :) = std_vec(j)*A_part_phi(j, :);
        D_thetaA(j, :) = std_vec_theta(j)*A(j, :);
        DA(j, :) = std_vec(j)*A(j, :);
    end
    c = zeros(M, 1);
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range) = sigma_j(index_j)*randn(Npix(index_j), 1);
        st = st+Npix(index_j);
    end
    neg_Y_theta_Gau_need(:, t) = (D_thetaA+DA_theta)*c*1e3;
    neg_Y_phi_Gau_need(:, t) = DA_phi*c*1e3;
    Y_Gau_need(:, t) = DA*c*1e3;
    
end

relative_energy_Gau_need = neg_Y_theta_Gau_need.^2+neg_Y_phi_Gau_need.^2;

save('sim_energy_Gau_need.mat', 'neg_Y_theta_Gau_need', 'neg_Y_phi_Gau_need',...
    'relative_energy_Gau_need', 'Y_Gau_need')