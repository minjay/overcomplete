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
burn_in = 1500;
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

samples = randsample(n_sample, T, true);

for t = 1:T
    
    i = samples(t);
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
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.2]);

cmax = 0;
for t = 1:9
    sim_whole_field = (large_scale+Y(:, t))/1e3;
    cmax = max(cmax, max(abs(sim_whole_field)));
end

phi_rot = phi+pi/2;
[x, y] = pol2cart(phi_rot, theta/pi*180);

for t = 1:9
    sim_whole_field = (large_scale+Y(:, t))/1e3;
    subplot(3, 3, t)
    cf = reshape(sim_whole_field, size(phi));
    vmag = linspace(min(cf(:)), max(cf(:)), 10);
    h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
    delete(h)
    shading flat
    caxis([-cmax cmax])
    text(-50, -50, sprintf('Min\n%2.1f',min(cf(:))),'FontName','times','Fontsize',10)
    text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf(:))),'FontName','times','Fontsize',10)
end
h = colorbar;
set(h, 'Position', [.85 .05 .05 .9]);