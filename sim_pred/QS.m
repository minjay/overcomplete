function value = QS(alpha, y, q)
% compute quantile score

value = ((y<q)-alpha)*(q-y);

end