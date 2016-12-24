function boxplot_group(data, g1, g2, name)
% data, g1, g2 are column vectors
% name is the name of the plot

bh = boxplot(data, {g1, g2}, 'colorgroup', g2, 'factorgap', 5, 'factorseparator', 1,...
    'Labels', {'nonGau', 'Gau', 'Matern','nonGau', 'Gau', 'Matern', 'nonGau', 'Gau', 'Matern'});
set(bh(6, :), 'LineWidth', 1.5)

set(gca, 'FontSize', 12)

annotation('textbox', [0.2 0.8 0.08 0.08], 'String', '\nu = 2.5', 'FitBoxToText', 'on',...
    'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
annotation('textbox', [0.475 0.8 0.08 0.08], 'String', '\nu = 3', 'FitBoxToText', 'on',...
    'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
annotation('textbox', [0.75 0.8 0.08 0.08], 'String', '\nu = 4', 'FitBoxToText', 'on',...
    'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
title(name)

end