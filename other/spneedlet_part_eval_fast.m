function psi_part = spneedlet_part_eval_fast(angle, B, j, bl_vector, P_prime, theta, phi,...
    xyz_xi, sqrt_lambda)

l_min = ceil(B^(j-1));
l_max = floor(B^(j+1));

n = length(theta);

psi_part = zeros(n, 1);
% dist_part: length n
if strcmp(angle, 'theta')
    dist_part = cos(theta).*cos(phi)*xyz_xi(1)+cos(theta).*sin(phi)*xyz_xi(2)...
    -sin(theta)*xyz_xi(3);
elseif strcmp(angle, 'phi')
    % here dist_part = correct dist_part/sin(theta)
    dist_part = -sin(phi)*xyz_xi(1)+cos(phi)*xyz_xi(2);
end

for i = 1:n
    psi_part(i) = sum(bl_vector(j+1, l_min:l_max).*P_prime(i, l_min+1:l_max+1));
end

psi_part = psi_part*sqrt_lambda.*dist_part;

end