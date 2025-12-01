addpath ../
clear
close all
rng(3);
figure('Position', [0, 0, 1100, 400])

% number of moments of the MGF that we compute
k = 80; 
% dimension of the basis for the approximation of the density
dim = 40; 

%% Linear fractional case
for c = [0.5, 0.7, 0.9]
    if c == 0.5
        color = 'm';
    elseif c == 0.7
        color = 'c';
    else
        color = 'k';
    end

    m = 1/(1-c);
    b = m * (1 - c)^2;
    p = [0; (b*c.^[0:k-2])'];
    p = p/sum(p);
    q = (1-b-c)/(c*(1-c));
    
    alpha = 1;
    true_sol = @(x) (1-q)^2*exp(-(1-q)*x);
    
    % Coefficients of phi(z)
    phi = newton_iter([1; -1], p, k, 100);
    if length(phi) < k
        phi = [phi; zeros(k-length(phi), 1)];
    end
    
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
    
    % Computation of the approximate density
    coeff = get_density_from_phi_a(phi, dim, alpha, beta, q);
    s = 1000;
    x = linspace(0, max(empirical_distr_W)*1.2, s);
    rho = evaluate_generalized_laguerre_fast(x, coeff, alpha, beta);
    
    % Plot
    subplot(1, 2, 1)
    plot(1:k, phi, strcat(color,'o'),'LineWidth', 2)
    hold on
    
    subplot(1, 2, 2) 
    tt = true_sol(x)';
    err = abs(rho-tt);
    semilogy(x, err, color, 'LineWidth', 2)
    hold on
end

subplot(1, 2, 1)
legend('c=0.5', 'c=0.7', 'c=0.9')
title('Coefficients of $\varphi$','Interpreter','latex');
set(gca,'fontsize', 16) 
subplot(1, 2, 2)
legend('c=0.5', 'c=0.7', 'c=0.9')
title('Error','Interpreter','latex')
set(gca,'fontsize', 16) 
print('-depsc', '-r600', 'testLF.eps')
hold off
