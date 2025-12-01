addpath ../
clear
close all
rng(0);

%% Example 6
p = [ones(10,1)/10];

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
T = 1000000; % Number of simulations
maxtime = 8; % Approximating W = W_maxtime

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
c = get_density_from_phi_a(phi(1:80), 40, alpha, beta, q);
s = 1000;
x = linspace(0, max(empirical_distr_W)*1.2, s);
rho = evaluate_generalized_laguerre_fast(x, c, alpha, beta);

%% Plot
hold on
plot(x, rho,'LineWidth', 2)
% plot(x, rho10,'LineWidth', 3)
%plot(x, rho5,'LineWidth', 3)

param = [1.8750699402245776, 1.538063418898302, 3.1295931802186137]; %20
a = param(1);
d = param(2);
p = param(3);
y = gengamma_pdf(x,a, d, p);

plot(x, y*(1-q), ':', 'LineWidth', 2)

param = [2.0362876789083364, 1.2338737165301643, 3.478181554698073]; % 10
a = param(1);
d = param(2);
p = param(3);
y = gengamma_pdf(x,a, d, p);

plot(x, y*(1-q), '--', 'LineWidth', 2)

param = [2.0893483732293117, 1.1744691568444119, 3.653307804549523]; % 5

a = param(1);
d = param(2);
p = param(3);
y = gengamma_pdf(x,a, d, p);

plot(x, y*(1-q), '-.', 'LineWidth', 2)

ylim([0, 0.75])
legend('Empirical pdf', 'Laguerre', 'Gamma (20)', 'Gamma (10)',...
    'Gamma (5)')
set(gca,'fontsize', 16) 
print('-depsc', '-r600', 'test6.eps')
hold off
