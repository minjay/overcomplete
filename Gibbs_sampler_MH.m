function post_samples = Gibbs_sampler_MH(A, Y, b_mat, fj_sq, nu, tau_sq_inv, V_inv, eta, tau_eta_sq, sigma_eta_sq, T, burn_in, thin, n_report)

[N, M] = size(A);
r = length(eta);
sample_size = floor((T-burn_in)/thin);
post_samples_c = zeros(M, sample_size);
post_samples_V_inv = zeros(M, sample_size);
post_samples_tau_sq_inv = zeros(1, sample_size);
post_samples_eta = zeros(r, sample_size);
flag = 0;
acc_times = 0;
for t = 1:T 
    
    if flag==0
        HatA = diag(exp(b_mat'*eta))*A;
    else
        HatA = HatA_star;
        flag = 0;
    end
    
    % sample c (slow)
    z = randn(M, 1);
    Sigma_inv = tau_sq_inv*(HatA'*HatA)+diag(V_inv);
    R = chol(Sigma_inv);
    z = z+R'\(HatA'*Y)*tau_sq_inv;
    c = R\z;
    
    % sample V
    shape = (nu+1)/2;
    scale = 2./(c.^2+nu*fj_sq);
    V_inv = gamrnd(shape, scale);
    
    % sample tau
    shape = N/2;
    HatAc = HatA*c;
    quad_form = (Y-HatAc)'*(Y-HatAc);
    scale = 2/quad_form;
    tau_sq_inv = gamrnd(shape, scale);
    
    % sample eta
    eta_star = mvnrnd(eta, sigma_eta_sq*eye(r))';
    f1 = tau_sq_inv*quad_form/2+eta'*eta/2/tau_eta_sq;
    Hat_star = diag(exp(b_mat'*eta_star));
    HatA_star = Hat_star*A;
    HatAc_star = HatA_star*c;
    quad_form_star = (Y-HatAc_star)'*(Y-HatAc_star);
    f2 = tau_sq_inv*quad_form_star/2+eta_star'*eta_star/2/tau_eta_sq;
    ratio = exp(f1-f2);
    u = rand;
    if ratio>=u
        eta = eta_star;
        flag = 1;
        acc_times = acc_times+1;
    end 
    
    if mod(t, n_report)==0
        disp(['Sampled: ', num2str(t), ' of ', num2str(T)])
        disp(['Overall Metrop. Acceptance rate: ', num2str(floor(acc_times/t*1e4)/100), '%'])
    end
    
    % save
    t_diff = t-burn_in;
    if t_diff>0 && mod(t_diff, thin)==0
        index = t_diff/thin;
        post_samples_c(:, index) = c;
        post_samples_V_inv(:, index) = V_inv;
        post_samples_tau_sq_inv(index) = tau_sq_inv;
        post_samples_eta(:, index) = eta;
    end
end

post_samples = struct('c', post_samples_c, 'V_inv', post_samples_V_inv,...
    'tau_sq_inv', post_samples_tau_sq_inv, 'eta', post_samples_eta);

end
    