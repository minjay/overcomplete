load('data_EOF_regr.mat')

phi = phi/pi*180;
theta = (pi/2-theta)/pi*180;

theta_vec = theta(:);
phi_vec = phi(:);
theta_min = min(theta_vec);

lonlim = [0 360];
latlim = [theta_min 90];

index = 1500*(1:9);
worldmap(latlim, lonlim)
for i = 1:9
    textm(theta_vec(index(i)), phi_vec(index(i)), num2str(i))
end

figure
for i = 1:9
    subplot(3, 3, i)
    [f, xi] = ksdensity(resid_all(:, index(i)));
    plot(xi, f)
    hold on
    ff = normpdf(xi, mean(resid_all(:, index(i))), std(resid_all(:, index(i))));
    plot(xi, ff, 'r')
end

figure
for i = 1:9
    h = subplot(3, 3, i);
    qqplot(resid_all(:, index(i)));
    delete(findall(h,'Type','text'))
    xlabel('Theoretical Quantiles')
    ylabel('Sample Quantiles')
end

    