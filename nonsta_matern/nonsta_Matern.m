function value = nonsta_Matern(std_s, std_t, r, nu, a)

if r>0
    value = 2^(1-nu)/gamma(nu)*(a*r)^nu*besselk(nu, a*r);
else
    value = 1;
end

value = value*std_s*std_t;

end