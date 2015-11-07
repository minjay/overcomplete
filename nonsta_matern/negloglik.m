function f_value = negloglik(beta_all, r, b_mat, samples)

n = length(samples);
beta = beta_all(1:end-1);
tau = beta_all(end);
cov_mat = get_cov(beta, r, b_mat)+eye(n)*tau^2;

f_value = ecmnobj(samples', zeros(n, 1), cov_mat);

end