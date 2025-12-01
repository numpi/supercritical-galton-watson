function y = evaluate_generalized_laguerre_fast(x, c, alpha, beta)
% y = evaluate_generalized_laguerre_fast(x, c, alpha, beta)
% Computes the sum of c_j (beta x)^(alpha-1) L_j(beta x) exp(-beta x)
% where L_j is the j-th generalized Laguerre polynomial with parameter
% alpha.

    a = alpha-1;
    x = x(:);
    s = length(x);

    Lminus2 = ones(s, 1);
    Lminus1 = ones(s, 1)*(1+a) - beta*x;
    y = c(1)*Lminus2 + c(2)*Lminus1;
    for i = 3:length(c)
        Lnew = ((2*i-3+a-beta*x).*Lminus1 - (i-2+a)*Lminus2)/(i-1);
        y = y + c(i)*Lnew;
        Lminus2 = Lminus1;
        Lminus1 = Lnew;
    end
    
    y = y .* exp(-beta*x) .* ((beta*x).^a);
end