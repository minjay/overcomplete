clear

% R is the radius of the ionosphere

load('theta_phi_R.mat')
load('ns.mat')
load('ns_deriv.mat')
load('mat_A.mat')
load('mat_A_part.mat')
load('post_samples_real_exp3.mat')

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
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(length(theta_vec), 1) b_mat];

m = size(b_mat, 2);

% posterior samples of std_vec
eta = beta_hat(1:m)';
std_vec = exp(b_mat*eta);

sigma_j = [1 sqrt(beta_hat(m+1:end-1))];

% subset A
[N, M] = size(A);
A = A(361:(N-360), :);
[N, M] = size(A);

b_mat_theta = kron(b_mat_deriv, ones(size(theta, 1), 1));
b_mat_theta = [zeros(length(theta_vec), 1) b_mat_theta];
std_vec_theta = exp(b_mat*eta).*(b_mat_theta*eta);

T = 1000;
E_theta_wo_coef_Gau_need = zeros(N, T);
E_phi_wo_coef_Gau_need = zeros(N, T);
DA_theta = zeros(N, M);
DA_phi = zeros(N, M);
D_thetaA = zeros(N, M);

Z_sim = randn(M, T);

for t = 1:T
    
    if mod(t, 10)==0
        t
    end
    
    for j = 1:N
        DA_theta(j, :) = std_vec(j)*A_part_theta(j, :);
        DA_phi(j, :) = std_vec(j)*A_part_phi(j, :);
        D_thetaA(j, :) = std_vec_theta(j)*A(j, :);
    end
    c = zeros(M, 1);
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range) = sigma_j(index_j)*Z_sim(range, t);
        st = st+Npix(index_j);
    end
    E_theta_wo_coef_Gau_need(:, t) = (D_thetaA+DA_theta)*c*1e3;
    E_phi_wo_coef_Gau_need(:, t) = DA_phi*c*1e3;
    
end

% the "4" comes from the stretching
% sin(\theta')/sin(\theta)
R = 6.5*1e6;
E_theta_Gau_need = -4*E_theta_wo_coef_Gau_need/R;
factor = sin(theta_vec*4)./sin(theta_vec);
E_phi_Gau_need = -repmat(factor, 1, T).*E_phi_wo_coef_Gau_need/R;

save('sim_energy_Gau_need.mat', 'E_theta_Gau_need', 'E_phi_Gau_need')