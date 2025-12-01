addpath ../

N = 100; % Number of moments

tol = 0;
maxit = 27;

% Generate starting point
phi0 =[1, -1];
% Generate an offspring distribution 
fprintf("Generating the offspring distribution...")
p = ones(1,10);
p = p/sum(p);
fprintf("...done!\n")
% End generation of the offspring distribution
% Fixed point iteration for the coefficients of phi(z)
fprintf("Computing coefficients of phi with fixed point...")
tic;
[phi_fx, err_fx] = simple_fixedpoint(phi0, p, N, maxit, tol);
cpu_time_fx = toc;
[phi_fx2, err_fx2] = simple_fixedpoint(phi0, p, N, maxit, tol, true);
fprintf("...done!\n")

% end fixed point
% Newton iteration for the coefficients of phi(z)
fprintf("Computing coefficients of phi with Newton...")
tic;
[phi_nw, err_nw] = newton_iter(phi0, p, N, maxit, tol); 
cpu_time_nw = toc;
[phi_nw2, err_nw2] = newton_iter(phi0, p, N, maxit, tol, true); 
fprintf("...done!\n")
% end newton

fprintf('---------Results---------\n')
fprintf('Fixed-point\n')
fprintf('Residual = %.4e, \t Iterations = %.1f,\t Cpu time = %.3f s\n',...
err_fx(end), length(err_fx), cpu_time_fx);
fprintf('Newton\n')
fprintf('Residual = %.4e, \t Iterations = %.1f,\t Cpu time = %.3f s\n',...
err_nw(end), length(err_nw), cpu_time_nw);
fprintf('------------------------------\n')

semilogy([1:length(err_fx)], err_fx,'b', [1:length(err_nw)], err_nw, 'r', [1:maxit], eps*ones(1,maxit), '--k');
figure
semilogy([1:length(err_fx2)], err_fx2,'b', [1:length(err_nw2)], err_nw2, 'r', [1:maxit], eps*ones(1,maxit), '--k'),
print('-depsc', '-r600', 'test_fixed_vs_newton.eps')


dlmwrite('test1_err1.dat', [[1:maxit]', err_fx, err_nw], '\t');
dlmwrite('test1_err2.dat', [[1:maxit]', err_fx2, err_nw2], '\t');




