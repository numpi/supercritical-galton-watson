dim_basis = 2:2:50;
cond_square = [];
cond_rectangular = [];
cond_square_scaled = [];
cond_rectangular_scaled = [];
for S = dim_basis
    S
    
    % Square case: # moments = # basis elements
    C = pascal(S+1, 1);
    C(:, 2:2:end) = -C(:, 2:2:end);
    cond_square = [cond_square, cond(vpa(C))];

    rescaling_factor = max(abs(C), [], 2);
    C = rescaling_factor.\C;
    cond_square_scaled = [cond_square_scaled, cond(vpa(C))];

    % Rectangular case: # moments = 2 * # basis elements
    C = pascal(2*S+2, 1);
    C = C(:, 1:S+1);
    C(:, 2:2:end) = -C(:, 2:2:end);
    cond_rectangular = [cond_rectangular, cond(vpa(C))];

    rescaling_factor = max(abs(C), [], 2);
    C = rescaling_factor.\C;
    cond_rectangular_scaled = [cond_rectangular_scaled, cond(vpa(C))];
end

semilogy(dim_basis, cond_square, '-o', 'linewidth', 2)
hold on
semilogy(dim_basis, cond_square_scaled, '-o', 'linewidth', 2)
semilogy(dim_basis, cond_rectangular, '-o', 'LineWidth', 2)
semilogy(dim_basis, cond_rectangular_scaled, '-o', 'LineWidth', 2)
legend('square', 'scaled square', 'rectangular',...
    'scaled rectangular', 'Location','best')
xlabel('S = basis dimension')
ylabel('condition number')
set(gca, 'fontsize', 16)
print('-depsc', '-r600', 'condition_numbers.eps')
hold off


