function value = CRPS_term2(mu, sigma_sq) 
% check paper "Probabilistic Forecasting and Comparative Model Assessment
% Based on Markov Chain Monte Carlo Output" Appendix B for details
% mu and sigma_sq are 1-D vectors

diff_mu = bsxfun(@minus, mu, mu');
sum_sigma_sq = bsxfun(@plus, sigma_sq, sigma_sq');

value = mean(mean(fun_A(diff_mu, sum_sigma_sq)))/2;

end
    
    