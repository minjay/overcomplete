addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

clear

[MSPE_2dot5, MAE_2dot5, CRPS_2dot5, QS_95_2dot5, QS_05_2dot5] = post_assess_eta0('2dot5');
[MSPE_3, MAE_3, CRPS_3, QS_95_3, QS_05_3] = post_assess_eta0('3');
[MSPE_4, MAE_4, CRPS_4, QS_95_4, QS_05_4] = post_assess_eta0('4');

MSPE_all = [MSPE_2dot5 MSPE_3 MSPE_4];
MAE_all = [MAE_2dot5 MAE_3 MAE_4];
CRPS_all = [CRPS_2dot5 CRPS_3 CRPS_4];
QS_95 = [QS_95_2dot5 QS_95_3 QS_95_4];
QS_05 = [QS_05_2dot5 QS_05_3 QS_05_4];

save('post_assess_nu.mat', 'MSPE_all', 'MAE_all', 'CRPS_all', 'QS_95', 'QS_05')
