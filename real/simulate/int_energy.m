load('sim_energy.mat')
load('sim_energy_Gau_need.mat')
load('sim_energy_Matern.mat')
load('WHI_quad_cond.mat')
load('theta_phi_R.mat')
load('energy_large_scale.mat')

P_cond = all_Cond_N{1}(:);
P_cond = P_cond(361:(end-360));

E_theta_large_scale = E_theta_large_scale(:);
E_phi_large_scale = E_phi_large_scale(:);

T = 5e3;
E_theta_sum = repmat(E_theta_large_scale, 1, T)+E_theta;
E_phi_sum = repmat(E_phi_large_scale, 1, T)+E_phi;
energy = (E_theta_sum.^2+E_phi_sum.^2).*repmat(P_cond, 1, T);
E_theta_sum_Gau_need = repmat(E_theta_large_scale, 1, T)+E_theta_Gau_need;
E_phi_sum_Gau_need = repmat(E_phi_large_scale, 1, T)+E_phi_Gau_need;
energy_Gau_need = (E_theta_sum_Gau_need.^2+E_phi_sum_Gau_need.^2).*repmat(P_cond, 1, T);
E_theta_sum_Matern = repmat(E_theta_large_scale, 1, T)+E_theta_Matern;
E_phi_sum_Matern = repmat(E_phi_large_scale, 1, T)+E_phi_Matern;
energy_Matern = (E_theta_sum_Matern.^2+E_phi_sum_Matern.^2).*repmat(P_cond, 1, T);
energy_large_scale = (E_theta_large_scale.^2+E_phi_large_scale.^2).*P_cond;

% compute the area of each latitudinal band
theta_one = theta(1, :);
n_theta = length(theta_one);
area_theta = zeros(1, n_theta);
for i = 1:n_theta
    if i==1
        theta_lower = 0;
    else
        theta_lower = (theta_one(i-1)+theta_one(i))/2;
    end
    if i==n_theta
        theta_upper = pi/4;
    else
        theta_upper = (theta_one(i)+theta_one(i+1))/2;
    end
    area_theta(i) = areaquad(90-theta_lower/pi*180, -180, 90-theta_upper/pi*180, 180);
end

% compute integrated energy for non-Gaussian and Gaussian needlet models
R = 6.5*1e6;
tot_area = 4*pi*R^2;
T = size(energy, 2);
int_energy_need = zeros(1, T);
int_energy_Gau_need = zeros(1, T);
int_energy_Matern = zeros(1, T);
for t = 1:T
    energy_one = reshape(energy(:, t), size(phi));
    int_energy_need(t) = sum(mean(energy_one, 1).*area_theta);
    energy_Gau_need_one = reshape(energy_Gau_need(:, t), size(phi));
    int_energy_Gau_need(t) = sum(mean(energy_Gau_need_one, 1).*area_theta);
    energy_Matern_one = reshape(energy_Matern(:, t), size(phi));
    int_energy_Matern(t) = sum(mean(energy_Matern_one, 1).*area_theta);
end

% compute integrated energy for large scale component
energy_large_scale_one = reshape(energy_large_scale, size(phi));
int_energy_large_scale = sum(mean(energy_large_scale_one, 1).*area_theta)*tot_area/1e9;

int_energy_need = int_energy_need*tot_area/1e9;
int_energy_Gau_need = int_energy_Gau_need*tot_area/1e9;
int_energy_Matern = int_energy_Matern*tot_area/1e9;

% set up colormap
map = brewermap(3, 'Set1');

figure
hold on

xs = linspace(0, max(int_energy_need), 200);
width = xs(2)-xs(1);
[nele, xs] = hist(int_energy_Matern, xs);
b1 = bar(xs-width/2, nele/T/width, 1, 'FaceColor', map(1, :));
b1.FaceAlpha = 0.5;
[f, xs] = ksdensity(int_energy_Matern, xs);
plot(xs, f, 'r', 'LineWidth', 2)
[nele, xs] = hist(int_energy_Gau_need, xs);
b2 = bar(xs-width/2, nele/T/width, 1, 'FaceColor', map(2, :));
b2.FaceAlpha = 0.5;
[f, xs] = ksdensity(int_energy_Gau_need, xs);
plot(xs, f, 'b', 'LineWidth', 2)
[nele, xs] = hist(int_energy_need, xs);
b3 = bar(xs-width/2, nele/T/width, 1, 'FaceColor', map(3, :));
b3.FaceAlpha = 0.5;
[f, xs] = ksdensity(int_energy_need, xs);
plot(xs, f, 'g', 'LineWidth', 2)
axis tight
y_range = get(gca, 'ylim');
h = plot([int_energy_large_scale int_energy_large_scale], y_range, 'k', 'LineWidth', 2);
lg = legend([b1 b2 b3 h], {'Gau-Matern', 'Gau-need', 'nonGau-need', 'large-scale'},...
    'location', 'northeast');
PatchInLegend = findobj(lg, 'type', 'patch');
set(PatchInLegend, 'facea', 0.5)
legend boxoff
set(gca, 'FontSize', 14)
xlim([0 max(int_energy_need)/2])
xlabel('Integrated Joule heating rate (GW)')
ylabel('Density')
