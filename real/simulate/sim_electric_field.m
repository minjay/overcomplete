% R is the radius of the ionosphere

load('theta_phi_R.mat')
load('deriv_B_spline.mat')
load('mat_A.mat')
load('mat_A_part.mat')
load('post_samples_real_new_resid_j24_pen_B_spline_cubic_2knots_aurora_n4000_sigma_j_sq_0.01_0.0002_nu4_long_init0.mat')

% set seed
rng(1)

% parameter setting
% d.f.
nu = 4;
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

% discard burn-in period
burn_in = 500;
post_samples_eta = post_samples.eta(:, burn_in+1:end);
post_samples_sigma_j = sqrt([1; 0.01; 0.0002]);
post_samples_tau = 1./sqrt(post_samples.tau_sq_inv(burn_in+1:end));

% num of posterior samples
n_sample = size(post_samples_eta, 2);

% posterior samples of std_vec
std_vec = exp(b_mat*post_samples_eta);

% subset A
[N, M] = size(A);
A = A(361:(N-360), :);
[N, M] = size(A);

b_mat_theta = kron(bS, ones(360, 1));
% the first column should be zero
b_mat_theta(:, 1) = 0;
std_vec_theta = exp(b_mat*post_samples_eta).*(b_mat_theta*post_samples_eta);

T = 5e3;
neg_Y_theta = zeros(N, T);
neg_Y_phi = zeros(N, T);
Y = zeros(N, T);
DA_theta = zeros(N, M);
DA_phi = zeros(N, M);
D_thetaA = zeros(N, M);
DA = zeros(N, M);

for t = 1:T
    
    i = n_sample;
    for j = 1:N
        DA_theta(j, :) = std_vec(j, i)*A_part_theta(j, :);
        DA_phi(j, :) = std_vec(j, i)*A_part_phi(j, :);
        D_thetaA(j, :) = std_vec_theta(j, i)*A(j, :);
        DA(j, :) = std_vec(j, i)*A(j, :);
    end
    c = zeros(M, 1);
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range) = post_samples_sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
        st = st+Npix(index_j);
    end
    neg_Y_theta(:, t) = (D_thetaA+DA_theta)*c*1e3;
    neg_Y_phi(:, t) = DA_phi*c*1e3;
    Y(:, t) = DA*c*1e3;
    
end

% the "4" comes from the stretching
relative_energy = (4*neg_Y_theta).^2+neg_Y_phi.^2;

save('sim_energy.mat', 'relative_energy', 'Y')