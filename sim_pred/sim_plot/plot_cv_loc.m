clear

name = '2dot5';

load(['data_sim_', name, '.mat'])

seed = 1;

rng(seed)

width = pi/9;

% sampling
[index, index_region] = rand_sampler(theta_vec, phi_vec, width);

N = length(theta_vec);
index_pred = setdiff(1:N, index);
index_pred_out_region = setdiff(index_pred, index_region);
index_pred_region = index_region;

scatter(phi_vec(index), theta_vec(index), 100, 'b.')
hold on
scatter(phi_vec(index_pred_out_region), theta_vec(index_pred_out_region), 100, 'g.')
scatter(phi_vec(index_pred_region), theta_vec(index_pred_region), 100, 'r.')
axis tight
axis equal

xlabel('Longitude (rad)')
ylabel('Co-latitude (rad)')
set(gca, 'FontSize', 12)
