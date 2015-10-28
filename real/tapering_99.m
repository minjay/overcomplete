% do tapering

load('data_99.mat')

resid_raw = resid;

theta_vec = theta(:)/pi*180;
n = length(theta_vec);

% 25 degree co-latitude
boundary = 25;

% std parameter
sigma_window = (max(theta_vec)-boundary)/1;

% get weights
index = find(theta_vec>boundary);
weights = ones(n, 1);
weights(index) = exp(-(theta_vec(index)-boundary).^2/2/sigma_window);

% plot weights
plot(theta_vec, weights, '.')

% tapering
resid = resid.*weights';

% plot
plot_pot_lite(reshape(resid, size(phi)), phi, theta, 1000, max(abs(resid)))

diff = resid_raw-resid;
plot_pot(reshape(diff, size(phi)), phi, theta, 1000)

std_vec_raw = std(reshape(resid_raw, size(phi)));
std_vec = std(reshape(resid, size(phi)));
plot(theta(1, :)/pi*180, std_vec_raw, '.')
hold on
plot(theta(1, :)/pi*180, std_vec, 'r.')

% save
save('data_99_taper.mat', 'resid', 'theta', 'phi')
