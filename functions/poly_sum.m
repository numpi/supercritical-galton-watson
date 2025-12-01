function h = poly_sum(f, g)
	d1 = length(f);
	d2 = length(g);
	h = zeros(max(d1, d2), 1);
	h(1:d1) = f(:);
	h(1:d2) = h(1:d2) + g(:);
end
