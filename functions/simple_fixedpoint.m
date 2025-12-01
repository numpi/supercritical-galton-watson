function [phi, err] = simple_fixedpoint(phi0, P, n, k, tol, do_check_poincare)
	if ~exist('tol', 'var')
		tol = 1e-14
	end
	P = P(:);
	m = sum([1:length(P)-1] .* P(2:end).');
	if ~exist('do_check_poincare', 'var')
		do_check_poincare = false;
	else
		do_check_poincare = true;	
	end
	if false && length(P) < n
		P = [P(:); zeros(n-length(P), 1)];
	end 
	if false && length(phi0) < n
		phi0 = [phi0(:); zeros(n-length(phi0), 1)];
	end 
	err = zeros(k, 1);
	d = m.^[0:n-1]';
	phi = phi0;
	for j=1:k
		phi_old = phi(:);
		phi = poly_comp_fft(P, phi, n);
		phi = d(1:length(phi)) .\ phi;
		phi(1) = 1;
		phi(2) = -1;
		
		if do_check_poincare
			err(j) = check_poincare(phi, P, n);
		else
			err(j) = norm(1-[phi_old(:); zeros(length(phi)-length(phi_old), 1)]./phi(:), inf);
		end	

		if err(j) < tol
			break
		end	
	end
	err = err(1:j);
end
