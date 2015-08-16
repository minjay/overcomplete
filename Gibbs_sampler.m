function post_samples = Gibbs_sampler(A, ATA, Y, ATY, fj_sq, nu, sigma_sq, tau_sq_inv, V_inv, T, burn_in, thin)

[N, M] = size(A);
sample_size = floor((T-burn_in)/thin);
post_samples_c = zeros(M, sample_size);
post_samples_V_inv = zeros(M, sample_size);
post_samples_sigma_sq = zeros(1, sample_size);
post_samples_tau_sq_inv = zeros(1, sample_size);
for t = 1:T 
    % sample c (slow)  
    z = randn(M, 1);
    Sigma_inv = tau_sq_inv*ATA+diag(V_inv);
    R = chol(Sigma_inv);
    z = z+R'\ATY*tau_sq_inv;
    c = R\z;
    
    % sample V
    shape = (nu+1)/2;
    scale = 2./(c.^2+nu*fj_sq*sigma_sq);
    V_inv = gamrnd(shape, scale);
    
    % sample sigma
    shape = nu*M/2;
    scale = 1/sum(fj_sq.*V_inv)*2/nu;
    sigma_sq = gamrnd(shape, scale);
    
    % sample tau
    shape = N/2;
    scale = 2/((Y-A*c)'*(Y-A*c));
    tau_sq_inv = gamrnd(shape, scale);
    
    % save
    t_diff = t-burn_in;
    if t_diff>0 && mod(t_diff, thin)==0
        index = t_diff/thin;
        post_samples_c(:, index) = c;
        post_samples_V_inv(:, index) = V_inv;
        post_samples_sigma_sq(index) = sigma_sq;
        post_samples_tau_sq_inv(index) = tau_sq_inv;
    end
end

post_samples = struct('c', post_samples_c, 'V_inv', post_samples_V_inv,...
    'sigma_sq', post_samples_sigma_sq, 'tau_sq_inv', post_samples_tau_sq_inv);

end
    