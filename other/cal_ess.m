function N = cal_ess(post_samples)
%CAL_ESS   Calculates the effective sample size based on the
%autocorrelation function (ACF).
%
%   N = cal_ess(post_samples);
%
% Inputs:
%   post_samples: posterior samples
% Outputs:
%   N: the effective sample size of the posterior samples
%
% Author: Minjie Fan, 2015

acf = autocorr(post_samples);

N = length(post_samples)/(1+2*sum(acf(2:end)));

end