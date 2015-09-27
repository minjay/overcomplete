function [pot_samples, theta_samples, phi_samples, index] = sampling_data(pot, theta, phi, n_samples, verbose)

theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec);

[pot_samples, index] = datasample(pot, n_samples, 'Replace', false,...
    'Weights', w);

theta_samples = theta_vec(index);
phi_samples = phi_vec(index); 

if verbose
    theta_vec = (pi/2-theta_vec)/pi*180;
    theta_min = min(theta_vec);
    lonlim = [0 360];
    latlim = [theta_min 90];

    worldmap(latlim, lonlim);
    lats = theta_vec(index);
    lons = phi_samples/pi*180;
    scatterm(lats, lons, [], pot_samples, 'filled');
    colorbar
    max_caxis = max(abs(pot_samples));
    caxis([-max_caxis max_caxis])
    
    setm(gca, 'ParallelLabel', 'off', 'MeridianLabel', 'off')
    textm(75, 300, '15\circ');
    textm(60, 300, '30\circ');
    textm(50, 177.5, '00')
    textm(52.5, 270, '06')
    textm(50, 2.5, '12')
    textm(47.5, 90, '18')
end

end
