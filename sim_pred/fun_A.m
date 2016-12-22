function value = fun_A(mu, sigma_sq)
% check paper "Probabilistic Forecasting and Comparative Model Assessment
% Based on Markov Chain Monte Carlo Output" Appendix B for details
% it supports vector/matrix inputs

sigma = sqrt(sigma_sq);
value = 2*sigma.*normpdf(mu./sigma)+mu.*(2*normcdf(mu./sigma)-1);

end