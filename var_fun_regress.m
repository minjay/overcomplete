n = 4;

mu = pi/4/(n+1)*(1:n);
rho = pi/4/(n+1)*2.5;

y = log_std_map_mean';
X = zeros(47, 1+n);
X(:, 1) = 1;
for i = 1:n
    X(:, i+1) = normpdf(theta, mu(i), rho/2);
end

beta = (X'*X)\(X'*y);

plot(90-theta/pi*180, exp(X*beta))
hold on
plot(90-theta/pi*180, exp(y), 'r')