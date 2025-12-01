%% Example of usage
% ----- Change the parameters in this section -----

% Set the offspring distrbution: it has to be a polynomial
p = [0.1; 0.3; 0.4; 0; 0.2];

% Set the number of moments of the moment generating
% function phi(z) of W to be computed
k = 80; 

% Set the dimension of the basis for the approximation
% of the density f(x) of the random variable W
dim = 40; 

% Set the maximum number of iteration in Newton's method
maxit = 100;

% Choose the parameters for the simulation (histogram)
T = 10000; % Number of simulations
maxtime = 15; % Approximating W = W_maxtime

%% The following line will compute the approximate density of W
[c, alpha, beta, q, x, y, phi, m, empirical_distr_W] = ...
    compute_density_W(p, k, dim, maxit, T, maxtime);

% Description of the outputs:
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

%% Plot the results
figure('Position',[0,0,500,300])
h = histogram(empirical_distr_W, 'Normalization','pdf');
h.FaceAlpha = 0.3;
hold on
plot(x, y, 'r','LineWidth', 3)
set(gca,'fontsize', 16) 
hold off