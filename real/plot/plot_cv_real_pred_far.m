clear

load('data_EOF_regr_new.mat')

resid = resid_all(1, :)'/1e3;

n = 1000;
theta_vec = theta(:);
phi_vec = phi(:);

width = pi/2;
lat_low = 20/180*pi;

rng(1)

% init weight vector w
w = sin(theta_vec*4);
% set the region of no data
w(theta_vec>=lat_low) = 0;
st = rand*2*pi;
en = st+pi/2;
% if part of the interval [st en] is outside of [0, 2*pi)
if en>=2*pi
    w(phi_vec>=st) = 0;
    w(phi_vec<=en-2*pi) = 0;
else
    w(phi_vec>=st & phi_vec<=en) = 0;
end

[pot_samples, index] = datasample(resid, n, 'Replace', false,...
    'Weights', w);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);

phi_samples_rot = phi_samples+pi/2;
[x_samples, y_samples] = pol2cart(phi_samples_rot, theta_samples/pi*180);

h = mypolar_grid([0 2*pi], [0 max(theta(:))/pi*180]);
delete(h)
hold on
scatter(x_samples(:), y_samples(:), 30, '.')
legend('Observed locations', 'Location', 'SouthOutside')
set(gca, 'FontSize', 12)
