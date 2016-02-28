function post_samples = Gibbs_sampler_AM3(model, data, params, tuning, options)
% GIBBS_SAMPLER_AM    The adaptive Metropolis-within-Gibbs sampler.
%
%   post_samples = Gibbs_sampler_AM(model, data, params, tuning, options);
%
% Inputs:
%   model - the model
%       model.A - the design matrix of the model, a N-by-M matrix
%       model.fj_sq - f(j)^2, a M-by-1 vector, the index corresponds to
%       psi_jk
%       model.b_mat - the design matrix of the function b, a N-by-(r+1)
%       matrix, each row corresponds to one location
%       model.nu - the number of degrees of freedom
%   data - the data
%       data.Y - the observed values, a N-by-1 vector
%       data.Npix - the number of spherical needlets at each frequency, a
%       (j_max-j_min+1)-by-1 vector
%   params - the initial values of the parameters
%       params.c - the initial value of c, a M-by-1 vector
%       params.V - the initial value of V^{-1}, a M-by-1 vector
%       params.eta = {eta_init, tau_sigma_sq, tau_eta_sq}
%           eta_init - the initial value of eta, a (r+1)-by-1 vector
%           tau_sigma_sq - the parameter tau_sigma_sq
%           tau_eta_sq - the parameter tau_eta_sq
%       params.tau - the initial value of 1/tau^2
%   tuning - the tuning parameters of the adaptive Metropolis
%       tuning.mu - the initial value of mu
%       (see Andrieu and Thoms, 2008, Algorithm 4)
%       tuning.Sigma - the initial value of Sigma
%       tuning.lambda - the initial value of lambda
%   options - the MCMC options
%       options.T - the number of MCMC iterations
%       options.burn_in - the length of the burn-in period
%       options.thin - the length of the thinning interval
%       options.n_report - the length of the interval to report progress
% Outputs:
%   post_samples - the posterior samples after discarding the samples in the
%   burn-in period and thinning
%       post_samples.c - the posterior samples of c, a M-by-sample_size
%       matrix
%       post_samples.V_inv - the posterior samples of V^{-1}, a
%       M-by-sample_size matrix
%       post_samples.tau_sq_inv - the posterior samples of 1/tau^2, a
%       1-by-sample_size vector
%       post_samples.eta - the posterior samples of eta, a
%       (r+1)-by-sample_size matrix
%
% Author: Minjie Fan, 2015

tic

% init
A = model.A;
b_mat = model.b_mat;
nu = model.nu;

Y = data.Y;
Npix = data.Npix;

c = params.c;
V_inv = params.V;
sigma_j_sq = [1; params.sigma_j_sq];
eta = params.eta;
eta_hat = params.eta_hat;
tau_sigma_sq = params.tau_sigma_sq;
tau_eta_sq = params.tau_eta_sq;
tau_sq_inv = params.tau;

mu = tuning.mu;
Sigma = tuning.Sigma;
lambda = tuning.lambda;

T = options.T;
burn_in = options.burn_in;
thin = options.thin;
n_report = options.n_report;

len_j = length(Npix);
st = zeros(len_j, 1);
en = zeros(len_j, 1);
for j = 1:len_j
    st(j) = sum(Npix(1:j))-Npix(j)+1;
    en(j) = sum(Npix(1:j));
end

[N, M] = size(A);

fj_sq = zeros(M, 1);
for j = 1:len_j
    range = st(j):en(j);
    fj_sq(range) = sigma_j_sq(j)*ones(Npix(j), 1);
end

r = length(eta)-1;

% optimal acceptance rate
% see Gelman et al., 1996
rates = [0.441 0.352 0.316 0.279 0.275 0.266 0.261 0.255 0.261 0.267 0.234];
if r+1<=10
    target_acc_rate = rates(r+1);
else
    target_acc_rate = rates(end);
end

sample_size = floor((T-burn_in)/thin);
post_samples_c = zeros(M, sample_size);
post_samples_V_inv = zeros(M, sample_size);
post_samples_sigma_j_sq = zeros(len_j, sample_size);
post_samples_tau_sq_inv = zeros(1, sample_size);
post_samples_eta = zeros(r+1, sample_size);

std_vec = exp(b_mat*eta);
DA = zeros(N, M);
for i = 1:N
    DA(i, :) = std_vec(i)*A(i, :);
end
acc_times = 0;

for t = 1:T 
    
    % sample c
    for j = 1:len_j  
        z = randn(Npix(j), 1);
        range = st(j):en(j);
        not_range = [1:st(j)-1, en(j)+1:M];
        DA_j = DA(:, range);
        DA_not_j = DA(:, not_range);
        Sigma_inv = tau_sq_inv*(DA_j'*DA_j)+diag(V_inv(range));
        R = chol(Sigma_inv);
        z = z+R'\(DA_j'*(Y-DA_not_j*c(not_range)))*tau_sq_inv;
        c(range) = R\z;
    end
   
    % sample V
    shape = (nu+1)/2;
    scale = 2./(c.^2+nu*fj_sq);
    V_inv = gamrnd(shape, scale);
    
    % sample sigma_j
    shape = nu*Npix(2:end)/2;
    scale = zeros(len_j-1, 1);
    for j = 1:len_j-1
        range = st(j+1):en(j+1);
        scale(j) = 1/sum(V_inv(range));
    end
    scale = 2/nu*scale;
    sigma_j_sq = [1; gamrnd(shape, scale)];
    for j = 1:len_j
        range = st(j):en(j);
        fj_sq(range) = sigma_j_sq(j)*ones(Npix(j), 1);
    end
    
    % sample tau
    shape = N/2;
    DAc = DA*c;
    quad_form = (Y-DAc)'*(Y-DAc);
    scale = 2/quad_form;
    tau_sq_inv = gamrnd(shape, scale);
    
    % sample eta
    eta_star = mvnrnd(eta, lambda*Sigma)';
    f1 = tau_sq_inv*quad_form/2+(eta(2:r+1)-eta_hat)'*(eta(2:r+1)-eta_hat)/2/tau_eta_sq+eta(1)^2/2/tau_sigma_sq;
    std_vec = exp(b_mat*eta_star);
    DA_star = zeros(N, M);
    for i = 1:N
        DA_star(i, :) = std_vec(i)*A(i, :);
    end
    DAc_star = DA_star*c;
    quad_form_star = (Y-DAc_star)'*(Y-DAc_star);
    f2 = tau_sq_inv*quad_form_star/2+(eta_star(2:r+1)-eta_hat)'*(eta_star(2:r+1)-eta_hat)/2/tau_eta_sq+eta_star(1)^2/2/tau_sigma_sq;
    ratio = exp(f1-f2);
    u = rand;
    % accept the new sample of eta
    if ratio>=u
        eta = eta_star;
        DA = DA_star;
        acc_times = acc_times+1;
    end 
    
    % adaptation step
    % gamma != 1/t to avoid being stuck at zero
    gamma = 1/(t+1);
    log_lambda = log(lambda)+gamma*(min([ratio 1])-target_acc_rate);
    lambda = exp(log_lambda);
    diff = eta-mu;
    mu = mu+gamma*diff;
    Sigma = Sigma+gamma*(diff*diff'-Sigma);
    
    % print to the screen
    if mod(t, n_report)==0
        if floor(t/n_report)==1
            disp('--------------------------------------------')
        end
        fprintf('Sampled: %d of %d\n', t, T)
        fprintf('Current Metropolis acceptance rate: %.2f%%\n',...
        acc_times/n_report*100)
        fprintf('Current lambda: %5f\n', lambda)
        disp('--------------------------------------------')
        acc_times = 0;
    end
    
    % save
    t_diff = t-burn_in;
    if t_diff>0 && mod(t_diff, thin)==0
        index = t_diff/thin;
        post_samples_c(:, index) = c;
        post_samples_V_inv(:, index) = V_inv;
        post_samples_sigma_j_sq(:, index) = sigma_j_sq;
        post_samples_tau_sq_inv(index) = tau_sq_inv;
        post_samples_eta(:, index) = eta;
    end
    
end

post_samples = struct('c', post_samples_c, 'V_inv', post_samples_V_inv,...
    'sigma_j_sq', post_samples_sigma_j_sq, 'tau_sq_inv',...
    post_samples_tau_sq_inv, 'eta', post_samples_eta);

toc

end
