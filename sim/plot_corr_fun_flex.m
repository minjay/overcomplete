function variance = plot_corr_fun_flex(B, sigma_j, j_min, j_max, nu, res, prop)

l_max = floor(B^(j_max+1));

dist_vec = linspace(-1, 1, res)';
Pl_mat = p_polynomial_value(res, l_max, dist_vec);
cov_vec = zeros(res, 1);
for j = j_min:j_max
    l_min = ceil(B^(j-1));
    l_max = floor(B^(j+1));
    tmp = zeros(res, 1);
    for l = l_min:l_max
        tmp = tmp+(fun_b(l/B^j, B))^2*(2*l+1)/(4*pi)*Pl_mat(:, l+1);
    end
    cov_vec = cov_vec+nu*(sigma_j(j-j_min+1))^2/(nu-2)*tmp;
end

plot(acos(dist_vec), cov_vec/cov_vec(dist_vec==1), prop, 'LineWidth', 2)
variance = cov_vec(dist_vec==1);

end