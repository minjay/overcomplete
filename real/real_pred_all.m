load('theta_phi.mat')
load('mat_A.mat')
% the output of real_fit.m
load('post_samples_real.mat')

% specify parameters
B = 2;
j_min = 2;
j_max = 4;
nu = 100;
alpha = 4;

theta_vec = theta(:);
phi_vec = phi(:);

[N, M] = size(A);

% non-stationary variance funcion
r = 4;
mu = pi/(r+1)*(1:r);
lambda = pi/(r+1)*2.5/2;
b_mat = zeros(N, r+1);
b_mat(:, 1) = 1;
for i = 2:r+1
    b_mat(:, i) = exp(-(theta_vec*4-mu(i-1)).^2/2/lambda^2);
end

% predict
T = 1000;
Ac = A*post_samples.c;
DAc = zeros(N, T);
for t = 1:T
    if mod(t, 10)==0
        t
    end
    eta = post_samples.eta(:, t);
    tau = 1/sqrt(post_samples.tau_sq_inv(t));
    std_vec = exp(b_mat*eta);
    for i = 1:N
        DAc(i, :) = std_vec(i)*Ac(i, :);
    end
end

% check posterior distribution
[f, xi] = ksdensity(DAc1(3000, :));
plot(xi, f)
hold on
[f, xi] = ksdensity(DAc2(3000, :));
plot(xi, f, 'r')

DAc_q1 = quantile(DAc1', [0.05, 0.5, 0.95]);
DAc_q2 = quantile(DAc2', [0.05, 0.5, 0.95]);
cmax = max(abs(DAc_q2(1,:)))*1e3;
plot_pot_lite(reshape(DAc_q1(1, :)*1e3, size(phi)), phi, theta, 1000, cmax)
figure
plot_pot_lite(reshape(DAc_q2(1,:)*1e3, size(phi)), phi, theta, 1000, cmax)


plot_pot(reshape(DAc_q1(3,:)*1e3, size(phi)), phi, theta, 1000)
figure
cmax = max(abs(DAc_q1(2,:)-DAc_q2(2,:)))*1e3;
plot_pot_lite(reshape(DAc_q1(2,:)-DAc_q2(2,:), size(phi))*1e3, phi, theta, 1000, cmax)
figure
plot_pot(reshape(DAc_q2(3,:)*1e3, size(phi)), phi, theta, 1000)

c = mean(post_samples.c, 2);
eta = mean(post_samples.eta, 2);
tau_sq = mean(1./sqrt(post_samples.tau_sq_inv));
std_vec = exp(b_mat*eta);
DA = zeros(N, M);
    for i = 1:N
        DA(i, :) = std_vec(i)*A(i, :);
    end
    DAc = DA*c;

Y_pred = DAc*1e3;

% plot fitted field
cmax = max([max(abs(Y_pred_nu3)) max(abs(Y_pred_nu100))]);
figure
plot_pot(reshape(Y_pred, size(phi)), phi, theta, 1000);
figure
plot_pot_lite(reshape(Y_pred_nu_4, size(phi)), phi, theta, 1000, cmax);

Y_pred_diff = Y_pred_nu3-Y_pred_nu100;

figure
plot_pot_lite(reshape(Y_pred_diff, size(phi)), phi, theta, 1000, max(abs(Y_pred_diff)));


load('data.mat')
figure
plot_pot(reshape(resid, size(phi)), phi, theta, 1000)

rng(1)

% sampling
N = 1e3;
[pot_samples, theta_samples, phi_samples, index] = sampling_data(resid,...
    theta, phi, N, 0);
% plot error field
Y_err = resid'-Y_pred_nu100;
figure
plot_pot_lite(reshape(Y_err, size(phi)), phi, theta, 1000, max(abs(Y_err)))

plot_pot(reshape(Y_pred_nu_3-Y_pred_nu_4, size(phi)), phi, theta, 1000)