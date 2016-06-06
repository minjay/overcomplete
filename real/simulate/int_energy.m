load('energy_both.mat')
load('theta_phi_R.mat')
load('energy_large_scale.mat')

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
int_energy = zeros(1, T);
int_energy_Gau_need = zeros(1, T);
for t = 1:T
    energy_one = reshape(energy(:, t)+energy_large_scale, size(phi));
    int_energy(t) = sum(mean(energy_one, 1).*area_theta);
    energy_Gau_need_one = reshape(energy_Gau_need(:, t)+energy_large_scale, size(phi));
    int_energy_Gau_need(t) = sum(mean(energy_Gau_need_one, 1).*area_theta);
end

% compute integrated energy for large scale component
energy_large_scale_one = reshape(energy_large_scale, size(phi));
int_energy_large_scale = sum(mean(energy_large_scale_one, 1).*area_theta)*tot_area/1e9;

int_energy = int_energy*tot_area/1e9;
int_energy_Gau_need = int_energy_Gau_need*tot_area/1e9;

% set up colormap
map = brewermap(2, 'Set1');

subplot(2, 1, 1)
xs = linspace(0, max(int_energy), 200);
histf(int_energy, xs,...
    'facecolor', map(1, :), 'facealpha', .5, 'edgecolor', 'none')
hold on
histf(int_energy_Gau_need, xs,...
    'facecolor', map(2, :), 'facealpha', .5, 'edgecolor', 'none')
axis tight
y_range = get(gca, 'ylim');
plot([int_energy_large_scale int_energy_large_scale], y_range, 'LineWidth', 2)
legalpha('nonGau-need','Gau-need', 'large-scale','location', 'northeast')
legend boxoff
set(gca, 'FontSize', 12)
xlim([0 max(int_energy)/2])
xlabel('Integrated Joule heating rate (GW)')

subplot(2, 1, 2)
int_energy_norm = (int_energy-mean(int_energy))/std(int_energy);
int_energy_Gau_need_norm = (int_energy_Gau_need-mean(int_energy_Gau_need))/...
    std(int_energy_Gau_need);
xs = linspace(min(int_energy_Gau_need_norm), max(int_energy_norm), 200);
[f, xi] = ksdensity(int_energy_norm, xs);
plot(xi, f, 'r', 'LineWidth', 2)
hold on
[f, xi] = ksdensity(int_energy_Gau_need_norm, xs);
plot(xi, f, 'b--', 'LineWidth', 2)
xlim([min(xs) max(xs)/2])
legend('nonGau-need','Gau-need','location', 'northeast')
legend boxoff
set(gca, 'FontSize', 12)

% save figure
% set white background
set(gcf, 'Color', 'w');
% 2.5 times resolution
export_fig int_energy_dist.png -m2.5
