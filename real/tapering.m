% do tapering

load('data.mat')

theta_vec = theta(:)/pi*180;
n = length(theta_vec);

% 25 degree co-latitude
boundary = 25;

% std parameter
sigma_window = (max(theta_vec)-boundary)/3;

% get weights
index = find(theta_vec>boundary);
weights = ones(n, 1);
weights(index) = exp(-(theta_vec(index)-boundary).^2/2/sigma_window);

% plot weights
plot(theta_vec, weights, '.')

% tapering
resid = resid'.*weights;

% plot
plot_pot_lite(reshape(resid, size(phi)), phi, theta, 1000, max(abs(resid)))

% save
save('data_taper.mat', 'resid', 'theta', 'phi')
