% load data
load('data_EOF_regr_new.mat')
resid = resid_all(1, :);

for i = 1:6
    subplot(3, 3, i)
    plot(post_samples.eta(i, :))
    axis tight
    title(['eta', num2str(i)])
end
subplot(3, 3, 7)
plot(post_samples.sigma_j_sq(2, :))
axis tight
hline = refline(0, 1e-2);
set(hline,'Color','r')
hline = refline(0, (1.5+2)/2*1e-2);
set(hline,'Color','r')
title('sigmajsq2')
subplot(3, 3, 8)
plot(post_samples.sigma_j_sq(3, :))
axis tight
hline = refline(0, 2e-4);
set(hline,'Color','r')
hline = refline(0, 3e-4);
set(hline,'Color','r')
title('sigmajsq3')
subplot(3, 3, 9)
plot(post_samples.tau_sq_inv)
axis tight
title('tausqinv')

load('mat_A.mat')
[N, M] = size(A);
A = A(361:N-360, :);
[N, M] = size(A);

theta_vec = theta(:);

% non-stationary variance function
knots = [0 0 0 0 40/180 80/180 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_vec*4);

b_mat(:, 1) = 1;

st = 1;
en = size(post_samples.eta, 2);
T = en-st+1;
Ac = A*post_samples.c(:, st:en);
DAc = zeros(N, T);
pred_err = zeros(T, 1);

% test prediction errors
for t = st:en
    std_vec = exp(b_mat*post_samples.eta(:, t));
    DAc(:, t-st+1) = std_vec.*Ac(:, t-st+1)*1e3;
    pred_err(t-st+1) = norm(DAc(:, t-st+1)-resid', 2);
end

% log plot
figure
semilogy(pred_err)

figure
hold on
std_vec = exp(b_mat*post_samples.eta);
% plot fitted std function
nu = 4;
for i = 1:100:size(post_samples.eta, 2)
    sigma_j = sqrt(post_samples.sigma_j_sq(:, i));
    variance = plot_corr_fun_flex(2, sigma_j, 2, 4, nu, 1000, 'r');
    std_vec(:, i) = std_vec(:, i)*sqrt(variance)*1000;
end
figure
plot(theta_vec, std_vec(:, 1:100:end))
