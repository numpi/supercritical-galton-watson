function [phi, err] = recursive_solution_fft(P, n, dP)
% Compute the solution of the Poincare equation by solving the non linear system via recursive substitution
	if isa(P, 'function_handle') && (~exist('dP', 'var') || ~isa(dP, 'function_handle'))
		error('If P is given as a function handle then you have to specify a function handle for dP as well')
	end
	if ~isa(P, 'function_handle')
		P = P(:);
		if ~exist('dP', 'var')
			dP = ([1:length(P)-1]') .* P(2:end);
		end	
	end			
	if isa(dP, 'function_handle')
		m = dP(1);
	else
		m = sum(dP);
	end
	phi = [1; -1; zeros(n-2, 1)];
	for j = 2:n-1
		if isa(P, 'function_handle')	
			B = toeplitz([phi(1:j); 0], [phi(1), zeros(1,j)]);
			Der = dP(B);
			Der = Der(end, end);
			Pol = P(B);
			Pol = Pol(end, 1);
			phi(j+1) = Pol/(m^j - Der);
		else
			Der = poly_comp_fft(dP, [phi(1:j); 0], j+1);
			Pol = poly_comp_fft(P, [phi(1:j); 0], j+1);
			phi(j+1) = Pol(end)/(m^j - Der(1));
		end	
    end	
    if ~isa(P, 'function_handle')
	    err = check_poincare(phi, P, n);
    else
        err = NaN;
    end
end
