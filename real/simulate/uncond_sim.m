clear

% load data
load('data_EOF_regr_new.mat')
% load pre-computed design matrix
load('mat_A.mat')
load('post_samples_real_reparam_nu3.mat')

% parameter setting
% d.f.
nu = 3;
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

% discard burn-in period
range = 2001:3000;
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

T = 4;
Y = zeros(N, T);
Y_comp = zeros(3, N, T);

rng(1)
samples = randsample(n_sample, T, true);

rng(2)
T_sim = trnd(nu, M, T);

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
        c(range) = post_samples_sigma_j(index_j, i)*T_sim(range, t);
        st = st+Npix(index_j);
        Y_comp(j-j_min+1, :, t) = DA(:, range)*c(range);     
    end
    Y(:, t) = (DA*c+post_samples_tau(i)*randn(N, 1));

end

phi_rot = phi+pi/2;
[x, y] = pol2cart(phi_rot, theta/pi*180);

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.2]);

cmax = max(max(abs(Y)));
% plot
for t = 1:T
    subplot(4, 4, t)
    cf = reshape(Y(:, t), size(phi));
    vmag = linspace(min(cf(:)), max(cf(:)), 10);
    h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
    delete(h)
    shading flat
    caxis([-cmax cmax])
    colormap(jet)
    title(num2str(t))
    text(-50, -50, sprintf('Min\n%2.1f',min(cf(:))),'FontName','times','Fontsize',10)
    text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf(:))),'FontName','times','Fontsize',10)
end
h = colorbar;
set(h, 'Position', [.85 1-0.2-0.04 .025 .2]);

for j = j_min:j_max
    cmax = max(max(abs(Y_comp(j-j_min+1, :, :))));
    for t = 1:T
        subplot(4, 4, (j-j_min+1)*T+t)
        cf = reshape(Y_comp(j-j_min+1, :, t), size(phi));
        vmag = linspace(min(cf(:)), max(cf(:)), 10);
        h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
        delete(h)
        shading flat
        caxis([-cmax cmax])
        colormap(jet)
        if j==j_min
            title('||')
        else
            title('+')
        end
        text(-50, -50, sprintf('Min\n%2.1f',min(cf(:))),'FontName','times','Fontsize',10)
        text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf(:))),'FontName','times','Fontsize',10)
    end
    h = colorbar;
    set(h, 'Position', [.85 1-(0.2+0.04)*(j-j_min+2) .025 .2]);
end

% save
save('uncond_sim.mat', 'Y', 'Y_comp')

print -painters -depsc uncond_sim_comp.eps