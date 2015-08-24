function [dist, psi] = get_psi_ss(B, j, k, theta, phi)
%GET_PSI   A wrapper of the function spneedlet_eval_fast. It evaluates the 
%spherical needlet with subscripts j and k at locations (theta, phi).
%
%   [dist, psi] = get_psi(B, j, k, theta, phi)
%
% Inputs:
%   B - the parameter
%   j, k - the subscripts
%   theta - the co-latitude of the locations, N-by-1 vector
%   phi - the longitude of the locations, N-by-1 vector
%
% Outputs:
%   dist - the inner products between the locations and the center of the
%   spherical needlet, N-by-1 vector
%   psi - the values of the spherical needlet at the locations, N-by-1
%   vector
%
% Author: Minjie Fan, 2015 

load('ss.mat')

j_max = j;
l_max = floor(B^(j_max+1));

bl_vector = get_bl_vector(B, j_max, l_max);

t = 2*floor(B^(j+1))+1;
grid_points = ss{degree_t==t};
sqrt_lambda = sqrt(4*pi/size(grid_points, 1));
x_xi = grid_points(k, 1);
y_xi = grid_points(k, 2);
z_xi = grid_points(k, 3);

[x, y, z] = sph2cart(phi, pi/2-theta, 1);

n = length(x);
dist = zeros(n, 1);

for i = 1:n
    dist(i) = sum([x(i) y(i) z(i)].*[x_xi y_xi z_xi]);
end

P = p_polynomial_value( n, l_max, dist );

psi = spneedlet_eval_fast(B, j, bl_vector, P, dist, sqrt_lambda);

end