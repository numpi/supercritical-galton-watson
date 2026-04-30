addpath ../
clear
close all
rng(1);
figure('Position', [0, 0, 600, 400])

% number of moments of the MGF that we compute
k = 80; 
% dimension of the basis for the approximation of the density
dim = 40; 

%% Linear fractional case
for c = [0.7, 0.9]
    m = 1/(1-c);
    b = m * (1 - c)^2;
    p = [0; (b*c.^[0:k-2])'];
    p = p/sum(p);
    q = (1-b-c)/(c*(1-c));
    
    alpha = 1;
    true_sol = @(x) (1-q)^2*exp(-(1-q)*x);
    
    % Coefficients of phi(z)
    phi_newton = newton_iter([1; -1], p, k, 100);
    if length(phi_newton) < k
        phi_newton = [phi_newton; zeros(k-length(phi_newton), 1)];
    end
    
    % Alternative: using fixed-point iteration
    % phi = simple_fixedpoint([1; -1], p, k, 1000);
    
    % Simulation and choice of beta
    T = 10000; % Number of simulations
    maxtime = 10; % Approximating W = W_maxtime
    if c == 0.9
        T = 1000;
        maxtime = 7;
    end
    
    empirical_distr_W = distr_W_from_simulation(p, T, maxtime);
    fig = figure(); set(fig,'visible','off');
    h = histogram(empirical_distr_W, 'Normalization','pdf');

    % Estimate of the beta parameter fron the histogram
    d = round(h.NumBins*0.3);
    data = [h.BinEdges(end-d:end)', h.Values(end-d:end)'];
    data(:, 2) = log(data(:, 2));
    ind_inf = find(data(:,2) == -inf);
    data(ind_inf, :) = [];
    beta = polyfit(data(:,1), data(:, 2), 1);
    beta = -beta(1);
    beta = max(1e-4, beta); % safeguard
    fprintf('Chosen beta: %1.2e\n', beta);
    close(fig)
    
    % Computation of the approximate density from Newton
    coeff_newton = get_density_from_phi_a(phi_newton, dim, alpha, beta, q);
    s = 1000;
    x = linspace(0, 15, s);
    rho_newton = evaluate_generalized_laguerre_fast(x, coeff_newton, alpha, beta);


    % Computation of coefficients from forward
    P = @(z) (1-c)*z / (eye(size(z))-c*z);
    dP = @(z) (1-c)*eye(size(z)) / ((eye(size(z))-c*z)^2);
    [phi_forward, ~] = recursive_solution_fft(P, k, dP);

    % Computation of the approximate density from forward
    coeff_forward = get_density_from_phi_a(phi_forward, dim, alpha, beta, q);
    s = 1000;
    x = linspace(0, 15, s);
    rho_forward = evaluate_generalized_laguerre_fast(x, coeff_forward, alpha, beta);
    
    % Plot
    if c == 0.7
        plot(1:k, phi_newton, 'sk','LineWidth', 1)
        hold on
        plot(1:k, phi_forward, 'om','LineWidth', 1)
    else
        plot(1:k, phi_newton, 'dc','LineWidth', 1)
        hold on
        plot(1:k, phi_forward, 'xb','LineWidth', 1)
    end
end

legend('c=0.7 Newton', 'c=0.7 forward', 'c=0.9 Newton', 'c=0.9 forward')
title('Coefficients of $\varphi$','Interpreter','latex');
set(gca,'fontsize', 16) 
print('-depsc', '-r600', 'testLF.eps')
hold off
