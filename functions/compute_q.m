function q = compute_q(p)
    if p(1) == 0
		q = 0;
    else
		p(2) = p(2)-1;
		r = roots(flip(p));
		epsilon = 1e-6;
		for i = 1:length(r)
		    if real(r(i))>=0 && real(r(i)) < 1-epsilon ...
		            && abs(imag(r(i))) < epsilon
		        q = real(r(i));
		    end
		end
    end
end