rng(1)

% parameter specification
nu = 3;
alpha = 4;

res = 200;

B = 2;
j_min = 2;
j_max = 4;
len_j = j_max-j_min+1;

% grid points
theta = linspace(0, pi, res/2);
phi = linspace(0, 2*pi, res);

[phi_mat, theta_mat] = meshgrid(phi, theta);
theta_vec = theta_mat(:);
phi_vec = phi_mat(:);

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_vec, phi_vec);
[N, M] = size(A); 

sigma_j = B.^(-alpha/2*(j_min:j_max));

% non-stationary variance funcion
r = 4;
mu = pi/(r+1)*(1:r);
lambda = pi/(r+1)*2.5/2;
b_mat = zeros(r+1, N);
b_mat(1, :) = 1;
for i = 2:r+1
    b_mat(i, :) = normpdf(theta_vec, mu(i-1), lambda);
end

eta = [1.5; randn(r, 1)];
A = diag(exp(b_mat'*eta))*A;

c = zeros(M, 1);
f_j = cell(len_j, 1);
st = 1;
for j = j_min:j_max
    index_j = j-j_min+1;
    range = st:st+Npix(index_j)-1;
    c(range) = sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
    f_j{index_j} = reshape(A(:, range)*c(range), res/2, res);
    st = st+Npix(index_j);
end

plot(c)

f = reshape(A*c, res/2, res);

[L, T] = meshgrid(phi-pi, pi/2-theta);
[HX, HY] = sph2hammer(L, T);

% plot
loc_x = [0 0.5 0 0.5];   
loc_y = [0.5 0.5 0 0];

figure
for j = j_min:j_max
    index_j = j-j_min+1;
    subplot('position', [loc_x(index_j) loc_y(index_j) 0.45 0.45])
    pcolor(HX, HY, f_j{index_j});
    shading interp
    axis equal
    axis tight
    axis off
    colorbar('southoutside')
    max_value = max(max(abs(f_j{index_j})));
    caxis([-max_value, max_value])
    title(['j = ', num2str(j)])
end

subplot('position', [loc_x(4) loc_y(4) 0.45 0.45])
pcolor(HX, HY, f);
shading interp
axis equal
axis tight
axis off
colorbar('southoutside')
max_value = max(max(abs(f)));
caxis([-max_value, max_value])
title('Sum')

% check non-Gaussianity
figure
qqplot(f(:))
