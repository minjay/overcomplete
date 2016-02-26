function index = rand_sampler_angle(theta_vec, phi_vec, angle)

index = find((phi_vec<=angle-pi/18 | phi_vec>=angle+pi/18) & theta_vec>=pi/6 & theta_vec<=5*pi/6);
n = 1e3;
index1 = randsample(index, n*0.9);
index = find((phi_vec<=angle-pi/18 | phi_vec>=angle+pi/18) & (theta_vec<pi/6 | theta_vec>5*pi/6));
index2 = randsample(index, n*0.1);
index = [index1; index2];

end
