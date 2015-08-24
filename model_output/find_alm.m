load('WHI_quad.mat')

all_Pot_N = all_Pot_N(1:720);
n = length(all_Pot_N);

% design matrix
X = zeros(n, numel(theta));
for i = 1:n
    X(i, :) = all_Pot_N{i}(:);
end

clear all_Pot_N

mean_map = mean(X, 1);
plot_pot(reshape(mean_map, size(phi)), phi, theta, 1000)

for i = 1:n
    X(i, :) = X(i, :)-mean_map;
end

[U, S, V] = svd(X, 'econ');

% find the 95% threshold
lambda = diag(S).^2;
var_exp = cumsum(lambda)./sum(lambda);
thres = find(var_exp>=0.99, 1, 'first');

% compute the residual fields
V_L = V(:, [1:thres]);
T_L = X*V_L;
X_L = T_L*V_L';
r = X-X_L;

bw = 180;
n = 4*bw^2;

% new theta and phi
phi_vec = 2*pi*(0:2*bw-1)/(2*bw);
theta_vec = pi*(2*(0:2*bw-1)+1)/4/bw;
[theta_interp, phi_interp] = meshgrid(theta_vec, phi_vec);

% project to the whole sphere
theta_global = theta*4;

% do interpolation
x_interp = interp2(theta_global, phi, reshape((X(1, :)), size(theta)),...
    theta_interp, phi_interp, 'spline', 0);

% crude comparison
subplot(1, 2, 1)
imagesc(reshape(X(1, :), size(phi)))
title('Before interpolation')
colorbar
subplot(1, 2, 2)
imagesc(x_interp)
title('After interpolation')
colorbar

x_interp_vec = x_interp(:);

filename = 'samples.dat';
fid = fopen(filename, 'w');

for i = 1:n
    fprintf(fid, '%.15f\n', x_interp_vec(i));
    fprintf(fid, '%.15f\n', 0);
end

alm = spharmonic_tran( filename, bw );

