function [Npix, grid_points, A] = get_A_ss(B, j_min, j_max, theta, phi)
%GET_A   Computes the design matrix A.
%
%   A = get_A(B, j_min, j_max, theta, phi, n_dist)
%
% Inputs:
%   B - the parameter
%   j_min - the minimal frequency
%   j_max - the maximal frequency
%   theta - the co-latitude of the locations, N-by-1 vector
%   phi - the longitude of the locations, N-by-1 vector
%   n_dist - the number of points on the fine grid
%
% Outputs:
%   A - the design matrix, N-by-M matrix
%
% Author: Minjie Fan, 2015

load('ss.mat')

N = length(theta);
[x, y, z] = sph2cart(phi, pi/2-theta, 1);
len_j = j_max-j_min+1;

Npix = zeros(len_j, 1);
grid_points = cell(len_j, 1);
for j = j_min:j_max
    index = j-j_min+1;
    t = 2*floor(B^(j+1))+1;
    grid_points{index} = ss{degree_t==t};
    Npix(index) = size(grid_points{index}, 1);
end

M = sum(Npix);
A = zeros(N, M);
dist = zeros(N, 1);
index_row = 0;

l_max = floor(B^(j_max+1));

bl_vector = get_bl_vector(B, j_max, l_max);

for j = j_min:j_max
    disp(['j = ', num2str(j), ' starts...'])
    index_j = j-j_min+1;
    l_max = floor(B^(j+1));
    sqrt_lambda = sqrt(4*pi/Npix(index_j));
    for k = 1:Npix(index_j)
        index_row = index_row+1;
        xyz_xi = grid_points{index_j}(k, :);
        for i = 1:N
            dist(i) = sum([x(i) y(i) z(i)].*xyz_xi);
            dist(i) = min(dist(i), 1);
            dist(i) = max(dist(i), -1);
        end
        P = p_polynomial_value( N, l_max, dist );
        A(:, index_row) = spneedlet_eval_fast(B, j, bl_vector, P, dist, sqrt_lambda);
    end
end

end