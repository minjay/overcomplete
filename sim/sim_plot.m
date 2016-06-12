subplot('position', [0.075 0.1 0.4 0.8])
plot_corr_fun(2, 2, 4, 3, 2.1, 10000, 'b-')
hold on
plot_corr_fun(2, 2, 4, 3, 3, 10000, 'r--')
plot_corr_fun(2, 2, 4, 3, 5, 10000, 'g-.')
legend({'\alpha=2.1', '\alpha=3', '\alpha=5'}, 'FontSize', 12)
xlabel('Great-circle distance (rad)')
ylabel('Correlation')
set(gca, 'FontSize', 12)
axis tight
axis square
title('(a)')

subplot('position', [0.575 0.1 0.4 0.8])
plot_corr_fun(2, 2, 4, 3, 3, 10000, 'b-')
hold on
plot_corr_fun(2, 3, 5, 3, 3, 10000, 'r--')
plot_corr_fun(2, 4, 6, 3, 3, 10000, 'g-.')
legend({'j=2-4', 'j=3-5', 'j=4-6'}, 'FontSize', 12)
xlabel('Great-circle distance (rad)')
ylabel('Correlation')
set(gca, 'FontSize', 12)
axis tight
axis square
title('(b)')