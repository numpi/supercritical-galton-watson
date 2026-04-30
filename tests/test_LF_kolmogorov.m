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
c = 0.7;
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
    
% Simulation and choice of beta
T = 10000; % Number of simulations 
maxtime = 10; % Approximating W = W_maxtime

empirical_distr_W = distr_W_from_simulation(p, T, maxtime);
fig = figure(); set(fig,'visible','off');
h = histogram(empirical_distr_W, 'Normalization','pdf');

% Estimate of the beta parameter fron the histogram
d = round(h.NumBins*0.4); 
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
s = 10000000;
x = linspace(0, 10, s);
rho_newton = evaluate_generalized_laguerre_fast(x, coeff_newton, alpha, beta);

% Computation of coefficients from forward
P = @(z) (1-c)*z / (eye(size(z))-c*z);
dP = @(z) (1-c)*eye(size(z)) / ((eye(size(z))-c*z)^2);
[phi_forward, ~] = recursive_solution_fft(P, k, dP);

% Computation of the approximate density from forward
coeff_forward = get_density_from_phi_a(phi_forward, dim, alpha, beta, q);
rho_forward = evaluate_generalized_laguerre_fast(x, coeff_forward, alpha, beta);

% Let's try: Kolmogorov distance between exact solution and simulations
k_dist_MC = zeros(1, maxtime);
T = 1000000; 
for gen = 1:maxtime
    fprintf("Simulating %d generations\n", gen);
    empirical_distr_gen = distr_W_from_simulation(p, T, gen);
    k_dist_MC(gen) = kolmogorov_dist_simulation(empirical_distr_gen);
end

% Kolmogorov distance between exact solution and our approximations
F_true = @(x) 1 - exp(-x);
F_ours_discretized = cumsum(rho_newton)/sum(rho_newton);
k_dist_ours = max(abs(F_true(x)' - F_ours_discretized));
F_forward_discretized = cumsum(rho_forward)/sum(rho_forward);
k_dist_forward = max(abs(F_true(x)' - F_forward_discretized));

semilogy(1:maxtime, k_dist_MC, '-*')
hold on
semilogy(1:maxtime, k_dist_ours * ones(1,maxtime), 'r', 'LineWidth', 2)
semilogy(1:maxtime, k_dist_forward * ones(1,maxtime), 'k--', 'LineWidth', 2)
xlabel('generation')
ylabel('Kolmogorov distance')
legend('Monte Carlo', 'Newton', 'forward')
set(gca,'fontsize', 16) 
print('-depsc', '-r600', 'testLF_distance0.7_1000000samples.eps')
hold off

function d = kolmogorov_dist_simulation(samples)
    samples = sort(samples);
    F_approx = @(x) sum(samples<x) / length(samples);
    F_true = @(x) 1-exp(-x);
    tt = linspace(0,samples(end)*2, 500);
    yy = zeros(length(tt), 1);
    for i = 1:length(yy)
        yy(i) = F_approx(tt(i));
    end
    d = max(abs(yy - F_true(tt)'));
end
