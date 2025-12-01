function err = check_poincare(phi, P, n, r, doplot)
	if ~exist('r','var')
		r = 1;
	end
	if ~exist('doplot','var')
		doplot = 0;
	end	
	P = P(:)';
	z = exp(2i*pi*[0:999]/1000)*r;
	m = sum([1:length(P)-1] .* P(2:end));

	temp = poly_comp_fft(P, phi, n);
	temp = temp(:)./(m.^[0:length(temp)-1]');
	err = 1 - [phi(:); zeros(length(temp)-length(phi), 1)]./temp;
	
	if doplot
		figure
		semilogy([1:1000], abs(err),'b', [1:1000], abs(polyval(phi(end:-1:1), m*z)), 'r', [1:1000], abs(polyval(P(end:-1:1), polyval(phi(end:-1:1), z))), 'c')
	end	
	err = norm(err, inf);
	
end
