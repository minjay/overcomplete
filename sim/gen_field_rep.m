rng(1)

% parameter specification
nu = 3;
alpha = 4;

res = 500;

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
[N, M] = size(A); 

sigma_j = B.^(-alpha/2*(j_min:j_max));

% the sample size needs to be large enough
% o.w., the normality test would be inaccurate
T = 1;
z = randn(M, T);
c = zeros(M, T);
for t = 1:T
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range, t) = sigma_j(index_j)*sqrt(nu)*z(range, t)./sqrt(chi2rnd(nu, Npix(index_j), 1));
        st = st+Npix(index_j);
    end
end

c_gauss = zeros(M, T);
for t = 1:T
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c_gauss(range, t) = sigma_j(index_j)*sqrt(nu/(nu-2))*z(range, t);
        st = st+Npix(index_j);
    end
end

st = 1;
[L, T] = meshgrid(phi-pi, pi/2-theta);
[HX, HY] = sph2hammer(L, T);
for j = j_min:j_max
    index_j = j-j_min+1;
    range = st:st+Npix(index_j)-1;
    f = A(:, range)*c(range, 1);
    f_gauss = A(:, range)*c_gauss(range, 1);
    st = st+Npix(index_j);
    
    % plot
    figure
    cmax = max([max(abs(f)) max(abs(f_gauss))]);
    subplot('position', [0.1 0.55 0.7 0.4])
    pcolor(HX, HY, reshape(f, size(phi_mat)));
    shading interp
    axis equal
    axis tight
    axis off
    caxis([-cmax cmax])
    title(['Simulated field, j=', num2str(j), ', nonGaussian, \nu=', num2str(nu)])
    subplot('position', [0.1 0.05 0.7 0.4])
    pcolor(HX, HY, reshape(f_gauss, size(phi_mat)));
    shading interp
    axis equal
    axis tight
    axis off
    caxis([-cmax cmax])
    title(['Simulated field, j=', num2str(j), ', Gaussian'])
    h = colorbar;
    set(h, 'Position', [0.85 0.05 0.05 0.9]);
end

f = A*c;
f_gauss = A*c_gauss;

figure
cmax = max([max(abs(f)) max(abs(f_gauss))]);
subplot('position', [0.1 0.55 0.7 0.4])
pcolor(HX, HY, reshape(f, size(phi_mat)));
shading interp
axis equal
axis tight
axis off
caxis([-cmax cmax])
title(['Simulated field, sum, nonGaussian', ', \nu=', num2str(nu)])
subplot('position', [0.1 0.05 0.7 0.4])
pcolor(HX, HY, reshape(f_gauss, size(phi_mat)));
shading interp
axis equal
axis tight
axis off
caxis([-cmax cmax])
title('Simulated field, sum, Gaussian')
h = colorbar;
set(h, 'Position', [0.85 0.05 0.05 0.9]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1)

T = 5000;
z = randn(M, T);
c = zeros(M, T);
for t = 1:T
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range, t) = sigma_j(index_j)*sqrt(nu)*z(range, t)./sqrt(chi2rnd(nu, Npix(index_j), 1));
        st = st+Npix(index_j);
    end
end

c_gauss = zeros(M, T);
for t = 1:T
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c_gauss(range, t) = sigma_j(index_j)*sqrt(nu/(nu-2))*z(range, t);
        st = st+Npix(index_j);
    end
end

f = A(1234, :)*c;
f_gauss = A(1234, :)*c_gauss;

subplot(1, 2, 1)
qqplot(f)
axis square
title(['Q-Q plot, nonGaussian', ', \nu=', num2str(nu)])
subplot(1, 2, 2)
qqplot(f_gauss)
axis square
title('Q-Q plot, Gaussian')
