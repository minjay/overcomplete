% plot sim_pred results

clear

[MAE_in_2dot5, MAE_out_2dot5, MSPE_in_2dot5, MSPE_out_2dot5, CRPS_in_2dot5, CRPS_out_2dot5,...
    QS_95_in_2dot5, QS_95_out_2dot5, QS_05_in_2dot5, QS_05_out_2dot5] = combine_results('2dot5');

[MAE_in_3, MAE_out_3, MSPE_in_3, MSPE_out_3, CRPS_in_3, CRPS_out_3,...
    QS_95_in_3, QS_95_out_3, QS_05_in_3, QS_05_out_3] = combine_results('3');

[MAE_in_4, MAE_out_4, MSPE_in_4, MSPE_out_4, CRPS_in_4, CRPS_out_4,...
    QS_95_in_4, QS_95_out_4, QS_05_in_4, QS_05_out_4] = combine_results('4');

MAE_in = [MAE_in_2dot5 MAE_in_3 MAE_in_4];
MAE_out = [MAE_out_2dot5 MAE_out_3 MAE_out_4];
MSPE_in = [MSPE_in_2dot5 MSPE_in_3 MSPE_in_4];
MSPE_out = [MSPE_out_2dot5 MSPE_out_3 MSPE_out_4];
CRPS_in = [CRPS_in_2dot5 CRPS_in_3 CRPS_in_4];
CRPS_out = [CRPS_out_2dot5 CRPS_out_3 CRPS_out_4];
QS_95_in = [QS_95_in_2dot5 QS_95_in_3 QS_95_in_4];
QS_95_out = [QS_95_out_2dot5 QS_95_out_3 QS_95_out_4];
QS_05_in = [QS_05_in_2dot5 QS_05_in_3 QS_05_in_4];
QS_05_out = [QS_05_out_2dot5 QS_05_out_3 QS_05_out_4];
MAE_in = MAE_in(:);
MAE_out = MAE_out(:);
MSPE_in = MSPE_in(:);
MSPE_out = MSPE_out(:);
CRPS_in = CRPS_in(:);
CRPS_out = CRPS_out(:);
QS_95_in = QS_95_in(:);
QS_95_out = QS_95_out(:);
QS_05_in = QS_05_in(:);
QS_05_out = QS_05_out(:);
g1 = [ones(60, 1); 2*ones(60, 1); 3*ones(60, 1)];
g2 = [ones(20, 3); 2*ones(20, 3); 3*ones(20, 3)];
g2 = g2(:);

% plot MAE_in
figure
boxplot_group(MAE_in, g1, g2, '(b) MAE (long-range)')

% plot MAE_out
figure
boxplot_group(MAE_out, g1, g2, '(a) MAE (short-range)')

% plot MSPE_in
figure
boxplot_group(MSPE_in, g1, g2, '(d) MSPE (long-range)')

% plot MSPE_out
figure
boxplot_group(MSPE_out, g1, g2, '(c) MSPE (short-range)')

% plot CRPS_in
figure
boxplot_group(CRPS_in, g1, g2, '(f) CRPS (long-range)')

% plot CRPS_out
figure
boxplot_group(CRPS_out, g1, g2, '(e) CRPS (short-range)')

% plot QS_05_in
figure
boxplot_group(QS_05_in, g1, g2, '(h) QS 5% (long-range)')

% plot QS_05_out
figure
boxplot_group(QS_05_out, g1, g2, '(g) QS 5% (short-range)')

% plot QS_95_in
figure
boxplot_group(QS_95_in, g1, g2, '(j) QS 95% (long-range)')

% plot QS_95_out
figure
boxplot_group(QS_95_out, g1, g2, '(i) QS 95% (short-range)')
