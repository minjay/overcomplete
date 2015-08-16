nu = 2;

alpha = 3;
sigma = 10;

res = 200;

theta = linspace(0, pi, res/2);
phi = linspace(0, 2*pi, res);

[phi_mat, theta_mat] = meshgrid(phi, theta);
theta_vec = theta_mat(:);
phi_vec = phi_mat(:);

B = 2;
j_min = 2;
j_max = 4;
n_dist = 1e3;

[Npix, A] = get_A_ss(B, j_min, j_max, theta_vec, phi_vec, n_dist);
[N, M] = size(A); 

sigma_j = sigma*B.^(-alpha*(2:4));

c = zeros(M, 1);
f_j = cell(j_max-j_min+1, 1);
index = 1;
for j = j_min:j_max
    range = index:index+Npix(j-j_min+1)-1;
    c(range) = sigma_j(j-j_min+1)*trnd(nu, Npix(j-j_min+1), 1);
    f_j{j-j_min+1} = reshape(A(:, range)*c(range), res/2, res);
    index = index+Npix(j-j_min+1);
end

plot(c)

f = reshape(A*c, res/2, res);

[L, T] = meshgrid(phi-pi, pi/2-theta);
[HX, HY] = sph2hammer(L, T);

figure
for j = 1:3
    subplot(2, 2, j)
    pcolor(HX, HY, f_j{j});
    shading flat
    axis equal
    axis tight
    axis off
    colorbar('southoutside')
    max_value = max(max(abs(f_j{j})));
    caxis([-max_value, max_value])
    title(['f', num2str(j)])
end

subplot(2, 2, 4)
pcolor(HX, HY, f);
shading flat
axis equal
axis tight
axis off
colorbar('southoutside')
max_value = max(max(abs(f)));
caxis([-max_value, max_value])
title('f1+f2+f3')