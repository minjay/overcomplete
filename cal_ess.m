function N = cal_ess(post_samples)
%CAL_ESS   Calculates the effective sample size based on the
%autocorrelation function (ACF).

acf = autocorr(post_samples);

N = length(post_samples)/(1+2*sum(acf(2:end)));

end