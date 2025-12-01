function c = get_density_from_phi_a(phi, dim, alpha, beta, q)
    % c = get_density_from_phi_a(phi, dim, alpha, beta, q)
    % Computes the coefficients of the density of W
    % INPUTS:
    % - phi: vector of coefficients of the Taylor expansion of the moment
    %   generating function of W
    % - dim: number of coefficients of the series expansion of the density
    %   of W that we want to compute
    % - alpha, beta: parameters in the series expansion of density of W
    % - q: probability of extinction

    a = alpha-1;
    momenti = abs(phi);
    momenti = momenti(:);
    momenti(1) = 1-q;
    num_momenti = length(phi);
    v = (1:num_momenti-1)./((1+a):(num_momenti+a-1));
    if abs(a) < 1e-8 
        % for numerical reasons to avoid division by almost zero
        momenti = momenti .* [1; cumprod(v)'];
    else
        momenti = (momenti ./ (a*gamma(a))) .* [1; cumprod(v)'];
    end
    momenti = momenti .* beta.^((1:num_momenti)');

    C = pascal(max(num_momenti, dim), 1);
    C = C(1:num_momenti, 1:dim);
    C(:, 2:2:end) = -C(:, 2:2:end);
    rescaling_factor = max(abs(C), [], 2);
    C = rescaling_factor.\C;
    b = momenti(1:num_momenti) ./ rescaling_factor; % RHS
    c = C \ b;
    c(2:2:end) = -c(2:2:end);
end
