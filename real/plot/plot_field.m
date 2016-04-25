load('data_EOF_regr_new.mat')
load('mat_A.mat')

theta_vec = theta(:);

A = A(:, 1:sum(Npix([1,2])));
N = size(A, 1);
A = A(361:(N-360), :);
[N, M] = size(A);

% non-stationary variance function
[~, ~, ~, X] = SCHA_regr(zeros(size(phi)), theta, phi, 3, 3, 4);
b_mat = X;

eta_mean = mean(post_samples.eta(:, 800:1000), 2);
std_vec = exp(b_mat*eta_mean);

for r = 1:5
    % predict
    T = 201;
    Ac = A*reshape(post_samples.c(:, r, 800:1000), M, T);
    DAc = zeros(N, T);
    for t = 1:T
        DAc(:, t) = std_vec.*Ac(:, t);
    end

    Y_pred = mean(DAc, 2)*1e3;

    figure
    plot_pot(reshape(Y_pred, size(phi)), phi, theta, 1000, max(abs(Y_pred)))
    
    Y_err = resid_all(1+3*(r-1), :)'-Y_pred;
    figure
    plot_pot_with_obs(reshape(Y_err, size(phi)), phi, theta, phi_samples, theta_samples, 1000)
end
