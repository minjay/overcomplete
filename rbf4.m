function f_value = rbf4(d, rho)

d = d/rho;
f_value = zeros(length(d), 1);
index = find(d<=1);
f_value(index) = (1-d(index)).^6.*(35*d(index).^2+18*d(index)+3)/3;

end