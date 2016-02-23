load('data_regr.mat')
load('resid_regr_all.mat')

phi = phi/pi*180;
theta = (pi/2-theta)/pi*180;

theta_vec = theta(:);
phi_vec = phi(:);
theta_min = min(theta_vec);

lonlim = [0 360];
latlim = [theta_min 90];

index = 1500*(1:9);
worldmap(latlim, lonlim)
scatterm(theta_vec(index), phi_vec(index), 'o', 'filled') 

for i = 1:9
    subplot(3, 3, i)
    [f, xi] = ksdensity(resid_all(:, index(i)));
    plot(xi, f)
    hold on
    ff = normpdf(xi, mean(resid_all(:, index(i))), std(resid_all(:, index(i))));
    plot(xi, ff, 'r')
end
    