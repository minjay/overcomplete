figure
for i = 1:25
    subplot(5, 5, i)
    ksdensity(relative_energy(i*500, :), 'support', 'positive')
end

figure
for i = 1:25
    subplot(5, 5, i)
    hist(relative_energy(i*500, :), 100)
end
