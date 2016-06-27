function index = rand_sampler(theta_vec, phi_vec, seed, width)

rng(seed)
lb = rand*(2*pi-width);
rb = lb+width;
index = find((phi_vec<=lb | phi_vec>=rb) & theta_vec>=pi/6 & theta_vec<=5*pi/6);
n = 1e3;
index1 = randsample(index, n*0.9);
index = find((phi_vec<=lb | phi_vec>=rb) & (theta_vec<pi/6 | theta_vec>5*pi/6));
index2 = randsample(index, n*0.1);
index = [index1; index2];

end
