% summarize PI and CP

clear
clc

[len_50_in_2dot5, len_50_out_2dot5, len_90_in_2dot5, len_90_out_2dot5,...
    cp_50_in_2dot5, cp_50_out_2dot5, cp_90_in_2dot5, cp_90_out_2dot5] = combine_results_PI('2dot5');
[len_50_in_3, len_50_out_3, len_90_in_3, len_90_out_3,...
    cp_50_in_3, cp_50_out_3, cp_90_in_3, cp_90_out_3] = combine_results_PI('3');
[len_50_in_4, len_50_out_4, len_90_in_4, len_90_out_4,...
    cp_50_in_4, cp_50_out_4, cp_90_in_4, cp_90_out_4] = combine_results_PI('4');

disp('50_out_2dot5')
round(mean(cp_50_out_2dot5, 1)*100, 1)
round(mean(len_50_out_2dot5, 1), 2)

disp('50_out_3')
round(mean(cp_50_out_3, 1)*100, 1)
round(mean(len_50_out_3, 1), 2)

disp('50_out_4')
round(mean(cp_50_out_4, 1)*100, 1)
round(mean(len_50_out_4, 1), 2)

disp('90_out_2dot5')
round(mean(cp_90_out_2dot5, 1)*100, 1)
round(mean(len_90_out_2dot5, 1), 2)

disp('90_out_3')
round(mean(cp_90_out_3, 1)*100, 1)
round(mean(len_90_out_3, 1), 2)

disp('90_out_4')
round(mean(cp_90_out_4, 1)*100, 1)
round(mean(len_90_out_4, 1), 2)

disp('50_in_2dot5')
round(mean(cp_50_in_2dot5, 1)*100, 1)
round(mean(len_50_in_2dot5, 1), 2)

disp('50_in_3')
round(mean(cp_50_in_3, 1)*100, 1)
round(mean(len_50_in_3, 1), 2)

disp('50_in_4')
round(mean(cp_50_in_4, 1)*100, 1)
round(mean(len_50_in_4, 1), 2)

disp('90_in_2dot5')
round(mean(cp_90_in_2dot5, 1)*100, 1)
round(mean(len_90_in_2dot5, 1), 2)

disp('90_in_3')
round(mean(cp_90_in_3, 1)*100, 1)
round(mean(len_90_in_3, 1), 2)

disp('90_in_4')
round(mean(cp_90_in_4, 1)*100, 1)
round(mean(len_90_in_4, 1), 2)
