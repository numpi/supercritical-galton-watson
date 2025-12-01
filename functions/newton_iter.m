function [phi, err] = newton_iter(phi0, P, n, k, tol, do_check_poincare)
	if ~exist('tol', 'var') 
		tol = 1e-14;
	end
	if ~exist('do_check_poincare', 'var') 
		do_check_poincare = false;
	else
		do_check_poincare = true;	
	end
	P = P(:);
	phi0 = phi0(:);
	m = sum([1:length(P)-1] .* P(2:end).');

	err = zeros(k, 1);
	d = m.^[0:n-1]';
	dP = [1:length(P)-1]' .* P(2:end);
	phi = phi0;
	for j=1:k
		phi_old = phi(:);
		col = poly_comp_fft(dP, phi, n);
        	thing = poly_sum(d(1:length(phi)) .* phi, -poly_comp_fft(P, phi, n));
        
        
		nn = max(length(col), length(thing));
		col = [col; zeros(nn -length(col), 1)];
		jac = diag(d(1:nn)) - toeplitz(col, [col(1), zeros(1, nn-1)]);

        % Let's update everything but the first two coefficients
        phismall = phi(3:end);
        jacsmall = jac(3:end, 3:end);
        thingsmall = thing(3:end);
	phismall = poly_sum(phismall, -jacsmall\thingsmall);
        phi = [1; -1; phismall];

		if do_check_poincare
			err(j) = check_poincare(phi, P, n); 
            fprintf("it = %d --> check_poincare = %1.2e\n", j, err(j));
		else
			err(j) = norm(1-[phi_old(:); zeros(length(phi)-length(phi_old), 1)]./phi(:), inf);
		end
			
		if err(j) < tol 
            fprintf("Number of iterations for Newton's method: %d\n", j);
			break
		end	
	end
	err = err(1:j);
end
