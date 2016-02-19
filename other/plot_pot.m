function plot_pot(Pot, phi, theta, res, max_caxis)

phi = phi/pi*180;
theta = (pi/2-theta)/pi*180;

theta_min = min(min(theta));

lonlim = [0 360];
latlim = [theta_min 90];

phi_interp = linspace(0, 360, res);
theta_interp = linspace(90, theta_min, res/4);

[theta_interp_mat, phi_interp_mat] = meshgrid(theta_interp, phi_interp);

Pot_interp = interp2(theta, phi, Pot, theta_interp_mat, phi_interp_mat, 'spline');

axes_pos = get(gca, 'Position');
axis off

ax1 = axes('Position', axes_pos, 'Visible', 'off');
x = axes_pos(1);
y = axes_pos(2);
width = axes_pos(3);
height = axes_pos(4);
axes('Position', [x+width/10, y+height/10, width-width/5, height-height/5]);

worldmap(latlim, lonlim);
pcolorm(theta_interp_mat, phi_interp_mat, Pot_interp/1e3);
h_c = colorbar;
ylabel(h_c, 'Potential Field [kV]');
max_caxis = max_caxis/1e3;
caxis([-max_caxis max_caxis])
setm(gca, 'ParallelLabel', 'off', 'MeridianLabel', 'off')
textm(75, 300, '15\circ');
textm(60, 300, '30\circ');
textm(50, 177.5, '00')
textm(52.5, 270, '06')
textm(50, 2.5, '12')
textm(47.5, 90, '18')

axes(ax1)
text(0.05, 0.1, ['min:', num2str(round(min(min(Pot_interp))/10)/1e2)])
text(0.05, 0.05, ['max:', num2str(round(max(max(Pot_interp))/10)/1e2)])

end
