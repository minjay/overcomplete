B = 2;
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
theta = zeros(N, 1);
phi = zeros(N, 1);
for i = 1:N
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end
[x, y, z] = sph2cart(phi, pi/2-theta, 1);

j_min = 2;
j_max = 4;
n_dist = 1e6;

[Npix, A] = get_A_ss(B, j_min, j_max, theta, phi, n_dist);
[N, M] = size(A);

mat = A'*A;
figure
imagesc(mat)
axis equal
axis tight
colorbar

mat2 = mat+diag(ones(M, 1));
tic
L = chol(mat2, 'lower');
toc

mat3 = diag(ones(M, 1));
tic
L2 = chol(mat3, 'lower');
toc


