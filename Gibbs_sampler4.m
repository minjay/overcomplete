function post_samples = Gibbs_sampler4(Npix, A, ATA, Y, ATY, fj_sq, nu, c, sigma_sq, tau_sq_inv, V_inv, T, burn_in, thin)

len_j = length(Npix);
st = zeros(len_j, 1);
en = zeros(len_j, 1);
for j = 1:len_j
    st(j) = sum(Npix(1:j))-Npix(j)+1;
    en(j) = sum(Npix(1:j));
end
    
[N, M] = size(A);
sample_size = floor((T-burn_in)/thin);
post_samples_c = zeros(M, sample_size);
post_samples_V_inv = zeros(M, sample_size);
post_samples_sigma_sq = zeros(1, sample_size);
post_samples_tau_sq_inv = zeros(1, sample_size);
for t = 1:T 
    if mod(t, 100)==0
        disp(['t = ', num2str(t), ' starts...'])
    end
    
    % sample c
    z = randn(M, 1);
    Sigma = 1./(tau_sq_inv*diag(ATA)+V_inv);
    diff_const = ATY-ATA*c;
    for k = 1:M
        diff_const = diff_const+ATA(:, k)*c(k);
        mu = tau_sq_inv*Sigma(k)*diff_const(k);    
        c(k) = mu+sqrt(Sigma(k))*z(k);
        diff_const = diff_const-ATA(:, k)*c(k);
    end
    
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
    