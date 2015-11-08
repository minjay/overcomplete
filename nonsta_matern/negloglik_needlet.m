function f_value = negloglik_needlet(beta_all, r, b_mat, samples, B, jj)

n = length(samples);
eta = beta_all(1:end-1);
tau = beta_all(end);
std_vec = exp(b_mat*eta);
corr_mat = zeros(n);
l_min = ceil(B^(jj-1));
l_max = floor(B^(jj+1));

index = 0;
r_vec = zeros((n+1)*n/2);
for j = 1:n
    for i = 1:j
        index = index+1;
        r_vec(index) = r(i, j);
    end
end
Pl_mat = p_polynomial_value(index, l_max, r_vec);

index = 0;
for j = 1:n
    j
    for i = 1:j
        index = index+1;
        value = 0;
        for l = l_min:l_max
            value = value+(fun_b(l/B^j, B))^2*(2*l+1)/(4*pi)*Pl_mat(index, l+1);
        end
        corr_mat(i, j) = value;
        corr_mat(j, i) = value;
    end
end

cov_mat = diag(std_vec)*corr_mat*diag(std_vec)+eye(n)*tau^2;

f_value = ecmnobj(samples', zeros(n, 1), cov_mat);

end