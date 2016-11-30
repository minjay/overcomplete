function [index, index_region] = rand_sampler(phi_vec, width)

lb = pi-width/2;
rb = pi+width/2;
index = find(phi_vec<=lb | phi_vec>=rb);
n = 500;
index = randsample(index, n);
index_region = find(phi_vec>lb & phi_vec<rb);

end
