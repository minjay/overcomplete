load('data_EOF_regr_new.mat')
resid = resid_all(1, :)';

figure

phi_rot = phi+pi/2;
[x, y] = pol2cart(phi_rot, theta/pi*180);

index = 3000*(1:5);

subplot('position', [0.05 0.35 0.6 0.6])
cmax = max(abs(resid/1e3));
cf = reshape(resid/1e3, size(phi));
vmag = linspace(min(cf(:)), max(cf(:)), 10);
h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
delete(h)
shading flat
caxis([-cmax cmax])
colormap(jet)
hold on
for i = 1:5
    text(x(index(i)), y(index(i)), num2str(i))
end

h = subplot('position', [0.65 1-0.3 0.2 0.2]);
qqplot(resid_all(:, index(1)));
delete(findall(h,'Type','text'))
ylabel('Sample Quantiles')
title('1')
axis square

h = subplot('position', [0.65 1-0.6 0.2 0.2]);
qqplot(resid_all(:, index(2)));
delete(findall(h,'Type','text'))
ylabel('Sample Quantiles')
title('2')
axis square

h = subplot('position', [0.65 1-0.9 0.2 0.2]);
qqplot(resid_all(:, index(5)));
delete(findall(h,'Type','text'))
xlabel('Theoretical Quantiles')
title('5')
axis square

h = subplot('position', [0.15 1-0.9 0.2 0.2]);
qqplot(resid_all(:, index(3)));
delete(findall(h,'Type','text'))
xlabel('Theoretical Quantiles')
ylabel('Sample Quantiles')
title('3')
axis square

h = subplot('position', [0.4 1-0.9 0.2 0.2]);
qqplot(resid_all(:, index(4)));
delete(findall(h,'Type','text'))
xlabel('Theoretical Quantiles')
title('4')
axis square
