addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

parpool(10)

name = '4';

T = 20;

% width of the region
% 24/12=2hr
width = pi/6;

parfor t = 1:T
    maxNumCompThreads(4);
    fit_sim_needlet(t, true, name, width);
end