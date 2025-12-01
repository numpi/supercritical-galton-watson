addpath ../
clear
close all
rng(0);

%% Example 5
p = [0; 0.3; 0.4; 0.2; 0.1];

% number of moments of the MGF that we compute
k = 80; 
% dimension of the basis for the approximation of the density
dim = 40; 

%% Coefficients of phi(z)
phi = newton_iter([1; -1], p, k, 100);
if length(phi) < k
    phi = [phi; zeros(k-length(phi), 1)];
end

%% The other parameters: q and alpha
q = compute_q(p);
mu = polyval(polyder(flip(p)), q);
m = polyval(polyder(flip(p)), 1);
fprintf("Value of m: %1.2e\n", m);
alpha = -log(mu) / log(m);
fprintf("Value of alpha: %1.2e\n", alpha);

%% Simulation and choice of beta
T = 100000; % Number of simulations
maxtime = 12; % Approximating W = W_maxtime

empirical_distr_W = distr_W_from_simulation(p, T, maxtime);

%%
figure('Position',[0, 0, 500, 300])
h = histogram(empirical_distr_W, 'Normalization','pdf');
h.FaceAlpha = 0.3;

% Estimate of the beta parameter fron the histogram
d = round(h.NumBins*0.3);
data = [h.BinEdges(end-d:end)', h.Values(end-d:end)'];
data(:, 2) = log(data(:, 2));
ind_inf = find(data(:,2) == -inf);
data(ind_inf, :) = [];
beta = polyfit(data(:,1), data(:, 2), 1);
beta = -beta(1);
fprintf('Chosen beta: %1.2e\n', beta);

%% Computation of the approximate density
c = get_density_from_phi_a(phi, 20, alpha, beta, q);
s = 1000;
x = linspace(0, max(empirical_distr_W)*1.2, s);
rho = evaluate_generalized_laguerre_fast(x, c, alpha, beta);

%% Plot
hold on
plot(x, rho,'LineWidth', 2)

% 20
param = [1.5529312744975565, 1.1887770591109867, 2.138336990025853];
a = param(1);
d = param(2);
p = param(3);
y = gengamma_pdf(x,a, d, p);

plot(x, y*(1-q), ':', 'LineWidth', 2)

% 10
param = [1.459414362701996, 1.3243823659287475, 2.045902913075773];
a = param(1);
d = param(2);
p = param(3);
y = gengamma_pdf(x,a, d, p);

plot(x, y*(1-q), '--', 'LineWidth', 2)

% 5
param = [1.2998579238969878, 1.5007965407838266, 1.869710834990371];

a = param(1);
d = param(2);
p = param(3);
y = gengamma_pdf(x,a, d, p);

plot(x, y*(1-q), '-.', 'LineWidth', 2)

ylim([-0.02, 1.02])
legend('Empirical pdf', 'Laguerre', 'Gamma (20)', 'Gamma (10)',...
    'Gamma (5)')
set(gca,'fontsize', 16) 
print('-depsc', '-r600', 'test5.eps')
hold off
