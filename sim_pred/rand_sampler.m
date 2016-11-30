function [index, index_region] = rand_sampler(phi, width)
% randomly split data into training and test data
% var index gives the index of training data
% var index_region gives the index of test data

lb = pi-width/2;
rb = pi+width/2;
index = find(phi<=lb | phi>=rb);
n = 500;
index = randsample(index, n);
index_region = find(phi>lb & phi<rb);

end
