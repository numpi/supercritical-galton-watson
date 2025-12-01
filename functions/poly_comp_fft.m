function h = poly_comp_fft(f, g, n)
	% compute the first nth coefficients of f(g(z)) mod z^n
	% f and g are given as f0 + f1 z + f2 z^2 + ...+ fs z^s 
	% Based on
	%
	% Kinoshita, Yasunori, and Baitian Li. "Power Series 
	% Composition in Near-Linear Time." 
	% arXiv preprint arXiv:2404.05177 (2024).
	
	f = f(:);
	m = length(f);
	
	% Shift the problem to match g(0) = 0
	if g(1) ~= 0
		f = poly_comp_dac(f, [g(1); 1]); % f(x) <-- f(x + g(0))
		g(1) = 0;						% g(x) <-- g(x) - g(0)
	end
	
	% Define auxiliary polynomials P and Q
	P = f(end:-1:1).';	 				% P(y) = y^m-1 f(y^-1)
	Q = [eye(length(g), 1), -g(:)]; 	% Q(x, y) = 1 - y g(x)
	
	d_max = (m-1)*(length(g)-1)+1; % degree of the composition
	% Call to the recursive part of the algorithm
	h = comp_fft_rec(P, Q, m-1, m, min(n, d_max));
end
%---------Recursive part of the procedure---------------------
function h = comp_fft_rec(P, Q, d, m, n)
	if n == 1
		h = P(d+1:min(m, length(P)));
	else
		% form sQ = Q(-x, y)
		sQ = Q;
		sQ(2:2:end, :) = -sQ(2:2:end, :); 
		
		% Compute V(x^2, y) = Q(x,y) * Q(-x, y) mod x^n mod y^m
		V = conv2(Q, sQ);
		V = V(1:min(n,size(V,1)), :);
		
		% Retrieve V(x,y) 
		V = V(1:2:end, :);
				
		% Recursive call
		e = max(0, d-size(Q, 2)+1);	
		W = comp_fft_rec(P, V, e, m, ceil(n/2));

		% Compute W(x^2, y) * Q(-x, y) mod x^n
		WW = zeros(2*size(W, 1)-1, size(W, 2));
		WW(1:2:end, :) = W;
		B = conv2(WW, sQ);
		B = B(1:min(n, size(B,1)), :);
		
		% Extract the coefficients from B
		h = B(:, d-e+1:m-e);
	end
end
