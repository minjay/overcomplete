function post_samples = Gibbs_sampler(A, ATA, Y, fj_sq, nu, sigma, tau, V, T, burn_in)

[N, M] = size(A);
sample_size = T-burn_in;
post_samples_c = zeros(M, sample_size);
post_samples_V = zeros(M, sample_size);
post_samples_sigma = zeros(1, sample_size);
post_samples_tau = zeros(1, sample_size);
for t = 1:T 
    t
    % sample c
    Sigma_inv = 1/tau^2*ATA+diag(1./V);
    mu = Sigma_inv\transpose(A)*Y/tau^2;
    R = chol(Sigma_inv);
    c = R\randn(M, 1)+mu;
    
    % sample V
    shape = (nu+1)/2;
    scale = 2./(c.^2+nu*fj_sq*sigma^2);
    V = 1./gamrnd(shape, scale);
    
    % sample sigma
    shape = nu*M/2;
    scale = 1/sum(fj_sq./V)*2/nu;
    sigma_sq = gamrnd(shape, scale);
    sigma = sqrt(sigma_sq);
    
    % sample tau
    shape = N/2;
    scale = 2/((Y-A*c)'*(Y-A*c));
    tau_sq = 1/gamrnd(shape, scale);
    tau = sqrt(tau_sq);
    
    % save
    if t>burn_in
        index = t-burn_in;
        post_samples_c(:, index) = c;
        post_samples_V(:, index) = V;
        post_samples_sigma(index) = sigma;
        post_samples_tau(index) = tau;
    end
end

post_samples = struct('c', post_samples_c, 'V', post_samples_V,...
    'sigma', post_samples_sigma, 'tau', post_samples_tau);

end
    