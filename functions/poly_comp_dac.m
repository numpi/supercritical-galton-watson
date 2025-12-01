function h = poly_comp_dac(f, g)
	% Composition of polynomials by means of divide-and-conquer 
	% efficient only when g has small degree: O(n log^2 n)  with deg(f) = d 
	f = f(:);
	g = g(:);
	d = length(f);
	if d == 1
		h = f(1);
	elseif d == 2
		h = f(2) * g(:);
		h(1) = h(1) + f(1);
	else	
		t = floor(d/2);
		% recursive calls
		h1 = poly_comp_dac(f(1:t), g);
		h2 = poly_comp_dac(f(t+1:end), g);
		% compute the tth power of g
		gtt = poly_power_fft(g, t);
		% compute the product gtt * h2
		gth2 = poly_mul_fft(gtt, h2);
		h = poly_sum(h1, gth2);
	end
end
