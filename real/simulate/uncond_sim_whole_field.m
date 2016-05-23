% load data
load('data_EOF_regr_new.mat')
% load pre-computed design matrix
load('mat_A.mat')
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

T = 9;
Y = zeros(N, T);
Y_comp = zeros(3, N, T);
rng(1)

for t = 1:T
    
    i = n_sample;
    DA = zeros(N, M);
    for j = 1:N
        DA(j, :) = std_vec(j, i)*A(j, :);
    end
    c = zeros(M, 1);
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range) = post_samples_sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
        st = st+Npix(index_j);
        Y_comp(j-j_min+1, :, t) = DA(:, range)*c(range)*1e3;     
    end
    Y(:, t) = (DA*c+post_samples_tau(i)*randn(N, 1))*1e3;

end

load('WHI_quad.mat')
load('theta_phi_R.mat')
whole_field = all_Pot_N{1}(:);
whole_field = whole_field(361:(end-360));
resid = resid_all(1, :)';
large_scale = whole_field-resid;

figure
plot_pot(reshape(large_scale, size(phi)), phi, theta, 1000, max(abs(large_scale)))

figure
for t = 1:9
    sim_whole_field = large_scale+Y(:, t);
    subplot(3, 3, t)
    plot_pot(reshape(sim_whole_field, size(phi)), phi, theta, 1000, max(abs(sim_whole_field)))
end
