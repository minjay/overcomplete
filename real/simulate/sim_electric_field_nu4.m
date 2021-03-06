clear

addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

% R is the radius of the ionosphere

load('theta_phi_R.mat')
load('ns.mat')
load('ns_deriv.mat')
load('mat_A.mat')
load('mat_A_part.mat')
load('post_samples_real_reparam_nu4.mat')

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
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(length(theta_vec), 1) b_mat];

range = 2001:3000;
% discard burn-in period
post_samples_eta = post_samples.eta(:, range);
post_samples_sigma_j = sqrt(post_samples.sigma_j_sq(:, range));
post_samples_tau = 1./sqrt(post_samples.tau_sq_inv(range));

% num of posterior samples
n_sample = size(post_samples_eta, 2);

% posterior samples of std_vec
std_vec = exp(b_mat*post_samples_eta);

% subset A
[N, M] = size(A);
A = A(361:(N-360), :);
[N, M] = size(A);

b_mat_theta = kron(b_mat_deriv, ones(size(theta, 1), 1));
b_mat_theta = [zeros(length(theta_vec), 1) b_mat_theta];
std_vec_theta = exp(b_mat*post_samples_eta).*(b_mat_theta*post_samples_eta);

T = 1000;
E_theta_wo_coef = zeros(N, T);
E_phi_wo_coef = zeros(N, T);
DA_theta = zeros(N, M);
DA_phi = zeros(N, M);
D_thetaA = zeros(N, M);

samples = randsample(n_sample, T, true);
T_sim = trnd(nu, M, T);

for t = 1:T
    
    if mod(t, 10)==0
        t
    end
    
    i = samples(t);
    for j = 1:N
        DA_theta(j, :) = std_vec(j, i)*A_part_theta(j, :);
        DA_phi(j, :) = std_vec(j, i)*A_part_phi(j, :);
        D_thetaA(j, :) = std_vec_theta(j, i)*A(j, :);
    end
    c = zeros(M, 1);
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range) = post_samples_sigma_j(index_j, i)*T_sim(range, t);
        st = st+Npix(index_j);
    end
    E_theta_wo_coef(:, t) = (D_thetaA+DA_theta)*c*1e3;
    E_phi_wo_coef(:, t) = DA_phi*c*1e3;
    
end

% the "4" comes from the stretching
% sin(\theta')/sin(\theta)
R = 6.5*1e6;
E_theta = -4*E_theta_wo_coef/R;
factor = sin(theta_vec*4)./sin(theta_vec);
E_phi = -repmat(factor, 1, T).*E_phi_wo_coef/R;

save('sim_energy_nu4.mat', 'E_theta', 'E_phi')
