function index = rand_sampler(phi_vec)

index = find(phi_vec<=pi-pi/18 | phi_vec>=pi+pi/18);
n = 1e3;
index = randsample(index, n);

end
