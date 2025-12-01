addpath ../
clear
close all
rng(0);
figure("Position",[0, 0, 600, 400])

%% Black robin example
pp = 0.1988;
b = zeros(7, 1);
for j = 0:5
    b(j+1) = nchoosek(5, j)*pp^j * (1-pp)^(5-j);
end
r = 1;
d = 1-0.6861;
p = zeros(7, 1);
p(1) = (1-r)*d + r*b(1)*d;
p(2) = (1-r)*(1-d) + r*(b(1)*(1-d) + b(2)*d);
for j = 3:7
    p(j) = r*(b(j-1)*(1-d) + b(j)*d);
end

p2 = [0.1036, 0.3551, 0.3448, 0.1553, 0.0366, 0.0044, 0.0002]';

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
maxtime = 15; % Approximating W = W_maxtime

empirical_distr_W = distr_W_from_simulation(p, T, maxtime);

%%
figure('Position', [0, 0, 500, 300])
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
c = get_density_from_phi_a(phi, dim, alpha, beta, q);
s = 1000;
x = linspace(0, max(empirical_distr_W)*1.2, s);
rho = evaluate_generalized_laguerre_fast(x, c, alpha, beta);

%% Plot
hold on
plot(x, rho, 'r','LineWidth', 3)

ylim([0, 1])
xlim([0, 6])
set(gca,'fontsize', 16) 
print('-depsc', '-r600', 'test_black_robin.eps')
hold off
