clear

name = '2dot5';

load(['data_sim_', name, '.mat'])
load(['post_samples_', name, '_1', '.mat'])

% convert theta and phi
theta = pi/2-theta;
% move the region to the center
phi = phi+pi;
% make sure phi is in [0, 2*pi)
phi(phi>=2*pi) = phi(phi>=2*pi)-2*pi;
phi(phi>pi) = phi(phi>pi)-2*pi;

[HX, HY] = sph2hammer(phi, theta);

N = length(theta);
index_pred = setdiff(1:N, index);
index_pred_out_region = setdiff(index_pred, index_region);
index_pred_region = index_region;

scatter(HX(index), HY(index), 20, 'bo')
hold on
scatter(HX(index_pred_out_region), HY(index_pred_out_region), 30, 'm+')
scatter(HX(index_pred_region), HY(index_pred_region), 30, 'r*')

% draw contour
th = linspace(-pi/2,pi/2,101);
lam = -pi+0*th;
[xh,yh] = sph2hammer(lam,th);
plot(xh,yh,'k');
lam = pi+0*th;
[xh,yh] = sph2hammer(lam,th);
plot(xh,yh,'k');

axis tight
axis off

legend('Observed locations', 'Short-range predicted locations', 'Long-range predicted locations',...
    'Location', 'SouthOutside')
set(gca, 'FontSize', 12)
