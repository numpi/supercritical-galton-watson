function [c, alpha, beta, q, x, y, phi, m, empirical_distr_W] = ...
    compute_density_W(p, k, dim, maxit, T, maxtime)
% INPUTS:
% - p is the coefficient vector of the polynomial that describes
%   the probability generating function of the offspring
%   distribution 
% - k is the number of coefficients of phi(z) that we wish to
%   compute. Default: k = 80
% - dim is the dimension of the basis that we wish to use for
%   the truncated series approximation of f(x), the density of W.
%   Default: dim = 40
% - maxit is the maximum number of iteration of Newton's method.
%   Default: maxit = 100
% - T is the number of simulations to be run in order to obtain an
%   estimate of the beta parameter. Default: T = 1000
% - maxtime is the time until the Galton-Watson processes are run
%   during the simulation. Default: maxtime = 12
%
% OUTPUTS:
% - c is the vector of coefficients in the series expansion of
%   the density function f(x) of W
% - alpha and beta are the parameters that appear in the 
%   series expansion of the density function f(x) of W
% - q is the extinction probability
% - y is the density function computed in the elements of
%   the vector x (equispaced between 0 and a suitable number)
% - phi is the vector that contains the first k coefficients
%   of the moment generating function of W (so a rescaled
%   version of the moments themselves)
% - m is the mean offsping of the distrbution (should be > 1
%   otherwise it is not a supercritical Galton-Watson process)
% - empirical_distr_W contains the data needed to plot the 
%   histogram of T simulations of the Galton-Watson process
%   until time "maxtime"

% Check inputs
if ~exist('k', 'var')
    k = 80;
end
if ~exist('dim', 'var')
    dim = 40;
end
if ~exist('maxit', 'var')
    maxit = 100;
end
if ~exist('T', 'var')
    T = 1000;
end
if ~exist('maxtime', 'var')
    maxtime = 12;
end

% Compute coefficients of phi(z)
phi = newton_iter([1; -1], p, k, maxit);
if length(phi) < k
    phi = [phi; zeros(k-length(phi), 1)];
end

% Computing q and alpha
q = compute_q(p);
mu = polyval(polyder(flip(p)), q);
m = polyval(polyder(flip(p)), 1);
alpha = -log(mu) / log(m);

% Simulation and choice of beta
empirical_distr_W = distr_W_from_simulation(p, T, maxtime);
h = histogram(empirical_distr_W, 'Normalization','pdf');
d = round(h.NumBins*0.3);
data = [h.BinEdges(end-d:end)', h.Values(end-d:end)'];
data(:, 2) = log(data(:, 2));
ind_inf = find(data(:,2) == -inf);
data(ind_inf, :) = [];
beta = polyfit(data(:,1), data(:, 2), 1);
beta = -beta(1);
close;

% Computation of the approximate density
c = get_density_from_phi_a(phi, dim, alpha, beta, q);
s = 1000;
x = linspace(0, max(empirical_distr_W)*1.2, s);
y = evaluate_generalized_laguerre_fast(x, c, alpha, beta);