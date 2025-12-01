function h = poly_power_fft(p, k)
	ss = k*(length(p)-1)+1;
	h = ifft(fft([p(:);zeros(ss-length(p), 1)]).^k);
end
