function index = rand_sampler_real(theta_vec)

index = find(theta_vec>=pi/6 & theta_vec<=5*pi/6);
n = 2000;
index1 = randsample(index, n*0.9);
index = find(theta_vec<pi/6 | theta_vec>5*pi/6);
index2 = randsample(index, n*0.1);
index = [index1; index2];

end
