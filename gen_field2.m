clear
load('std_map.mat')

nu = 3;

alpha = 3;
sigma = 10;

res = 200;

theta = linspace(0, pi, res/2);
phi = linspace(0, 2*pi, res);

[theta_mat, phi_mat] = meshgrid(theta, phi);
theta_vec = theta_mat(:);
phi_vec = phi_mat(:);

B = 2;
j_min = 2;
j_max = 4;
n_dist = 1e3;

[Npix, grid_points, A] = get_A_ss(B, j_min, j_max, theta_vec, phi_vec, n_dist);
[N, M] = size(A); 

c = zeros(M, 1);
f_j = cell(j_max-j_min+1, 1);
std_map_interp = interp2(theta_map*4, phi_map, reshape(std_map, size(theta_map)), theta_vec, phi_vec, 'spline');
index = 0;
for j = j_min:j_max
    range = sum(Npix(1:j-j_min+1))-Npix(j-j_min+1)+1:sum(Npix(1:j-j_min+1));
    for k = 1:Npix(j-j_min+1)
        index = index+1;
        xyz_xi = grid_points{j-j_min+1}(k, :);
        [phi_xi, theta_xi, ~] = cart2sph(xyz_xi(1), xyz_xi(2), xyz_xi(3));
        if phi_xi<0
            phi_xi = phi_xi+2*pi;
        end
        theta_xi = pi/2-theta_xi;
        sigma_jk = sigma*...
            B^(-alpha*j);
        c(index) = sigma_jk*trnd(nu);
    end
    f_j{j-j_min+1} = reshape(A(:, range)*c(range).*std_map_interp, res, res/2);
end

plot(c)

f = reshape(A*c.*std_map_interp, res, res/2);

for j = 1:3
    figure
    plot_pot(f_j{j}, phi_mat, theta_mat/4, 1000)
end

figure
plot_pot(f, phi_mat, theta_mat/4, 1000)
title('f1+f2+f3')