function h = poly_mul_fft(p, q)
	d1 = length(p);
	d2 = length(q);
	h = ifft(fft([p(:); zeros(d2-1,1)]) .* fft([q(:); zeros(d1-1,1)]));
end
