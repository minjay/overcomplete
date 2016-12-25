% load data
load('data_EOF_regr_new.mat')
% load pre-computed design matrix
load('mat_A.mat')
load('post_samples_real_exp3.mat')
load('post_samples_exp2.mat')

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
load('ns.mat')
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(length(theta_vec), 1) b_mat];

m = size(b_mat, 2)-1;
beta_z = beta_hat(1:end-1);
tau_z = beta_hat(end);
sigma_sq_z = [1 beta_z(m+2:end)];

% discard burn-in period
post_samples_eta = post_samples.eta(:, 2001:3000);
post_samples_sigma_j = sqrt(post_samples.sigma_j_sq(:, 2001:3000));
post_samples_tau = 1./sqrt(post_samples.tau_sq_inv(2001:3000));

% num of posterior samples
n_sample = size(post_samples_eta, 2);

% posterior samples of std_vec
std_vec = exp(b_mat*post_samples_eta);

% subset A
[N, M] = size(A);
A = A(361:(N-360), :);
[N, M] = size(A);

T = 9;
Y_z = zeros(N, T);
Y_t = zeros(N, T);
rng(1)

samples = randsample(n_sample, T, true);

rng(2)
Z_sim = randn(M, T);

rng(2)
T_sim = trnd(nu, M, T);

err_sim = randn(N, T);

for t = 1:T
    
    i = samples(t);
    DA = zeros(N, M);
    for j = 1:N
        DA(j, :) = std_vec(j, i)*A(j, :);
    end
    c_z = zeros(M, 1);
    c_t = zeros(M, 1);
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c_z(range) = sqrt(sigma_sq_z(index_j))*Z_sim(range, t);
        c_t(range) = post_samples_sigma_j(index_j, i)*T_sim(range, t);
        st = st+Npix(index_j);
    end
    Y_z(:, t) = (DA*c_z+tau_z*err_sim(:, t))*1e3;
    Y_t(:, t) = (DA*c_t+post_samples_tau(i)*err_sim(:, t))*1e3;

end

Y_sim_Gau_need = Y_z';
Y_sim_need = Y_t';

save('Y_sim_need.mat', 'Y_sim_Gau_need', 'Y_sim_need')

load('WHI_quad.mat')
load('theta_phi_R.mat')
whole_field = all_Pot_N{1}(:);
whole_field = whole_field(361:(end-360));
resid = resid_all(1, :)';
large_scale = whole_field-resid;

figure
plot_pot(reshape(large_scale, size(phi)), phi, theta, 1000, max(abs(large_scale)))

figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.2]);

cmax = 0;
for t = 1:9
    sim_whole_field = (large_scale+Y_t(:, t))/1e3;
    cmax = max(cmax, max(abs(sim_whole_field)));
end

phi_rot = phi+pi/2;
[x, y] = pol2cart(phi_rot, theta/pi*180);

for t = 1:9
    sim_whole_field = (large_scale+Y_t(:, t))/1e3;
    subplot(3, 3, t)
    cf = reshape(sim_whole_field, size(phi));
    vmag = linspace(min(cf(:)), max(cf(:)), 10);
    h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
    delete(h)
    shading flat
    caxis([-cmax cmax])
    colormap(jet)
    text(-50, -50, sprintf('Min\n%2.1f',min(cf(:))),'FontName','times','Fontsize',10)
    text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf(:))),'FontName','times','Fontsize',10)
end
h = colorbar;
set(h, 'Position', [.85 .05 .05 .9]);