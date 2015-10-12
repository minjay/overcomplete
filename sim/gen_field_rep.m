rng(1)

% parameter specification
nu = 3;
alpha = 4;

res = 100;

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

T = 1000;
c = zeros(M, T);
for t = 1:T
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range, t) = sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
        st = st+Npix(index_j);
    end
end
f = A*c;
