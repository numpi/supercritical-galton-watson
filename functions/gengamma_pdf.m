function y = gengamma_pdf(x, a, d, p)
    % Generalized Gamma PDF
    y = p / (a^d * gamma(d / p)) * x.^(d - 1) .* exp(-(x / a).^p);
end
