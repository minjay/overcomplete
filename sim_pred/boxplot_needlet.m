load('sim_pred_err_Matern_rep_2dot5.mat')
MSPE_out_2dot5 = MSPE_Matern_out_region_all;
MSPE_2dot5 = MSPE_Matern_region_all;
MAE_out_2dot5 = MAE_Matern_out_region_all;
MAE_2dot5 = MAE_Matern_region_all;

load('sim_pred_err_Matern_rep_3.mat')
MSPE_out_3 = MSPE_Matern_out_region_all;
MSPE_3 = MSPE_Matern_region_all;
MAE_out_3 = MAE_Matern_out_region_all;
MAE_3 = MAE_Matern_region_all;

load('sim_pred_err_Matern_rep_4.mat')
MSPE_out_4 = MSPE_Matern_out_region_all;
MSPE_4 = MSPE_Matern_region_all;
MAE_out_4 = MAE_Matern_out_region_all;
MAE_4 = MAE_Matern_region_all;

subplot(2, 2, 1)
boxplot([MSPE_out_2dot5 MSPE_out_3 MSPE_out_4], 'Labels', {'nu=2.5', 'nu=3', 'nu=4'})
ylabel('MSPE')
title('Out of region')
subplot(2, 2, 2)
boxplot([MSPE_2dot5 MSPE_3 MSPE_4], 'Labels', {'nu=2.5', 'nu=3', 'nu=4'})
ylabel('MSPE')
title('In region')
subplot(2, 2, 3)
boxplot([MAE_out_2dot5 MAE_out_3 MAE_out_4], 'Labels', {'nu=2.5', 'nu=3', 'nu=4'})
ylabel('MAE')
title('Out of region')
subplot(2, 2, 4)
boxplot([MAE_2dot5 MAE_3 MAE_4], 'Labels', {'nu=2.5', 'nu=3', 'nu=4'})
ylabel('MAE')
title('In region')
