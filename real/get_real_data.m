load('svd.mat')
load('theta_phi.mat')

t = 1;
resid = double(r(t, :));
save('data.mat', 'resid', 'theta', 'phi')