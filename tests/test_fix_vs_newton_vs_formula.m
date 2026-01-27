addpath ../

degree = [5, 10, 15, 20, 50];
avg_offsp = [1.1 1.25 1.5 2 3];
N = 100; % Number of moments
tol = 1e-8;
maxit = 300;
ntests = 10;
do_check_poincare = false;
% Data structure: for any pair (d,m) store the residual, number of iterations, and cpu time
data_fixp = zeros(length(degree), length(avg_offsp), 3);
data_newt = zeros(length(degree), length(avg_offsp), 3);
data_recu = zeros(length(degree), length(avg_offsp), 2);
i1=1;
for d = degree
	i2=1;
	for m = avg_offsp
		for j = 1:ntests
			% Generate starting point
			%phi0 = [1, -1, rand(1, N-2)-0.5];
			phi0 =[1, -1];
			% Generate an offspring distribution with p'(1) = m and with p poly of degree d
			fprintf("Generating the offspring distribution...")
			p = 0;
			while ([1:length(p)] .* p.') < m
				p = rand(1, d);
				p = p/sum(p);
			end
			p = p * m / ([1:length(p)] * p.');
			p = [1 - sum(p), p];
			fprintf("...done!\n")
			% End generation of the offspring distribution
			% Fixed point iteration for the coefficients of phi(z)
			fprintf("Computing coefficients of phi with fixed point...")
			tic;
			[phi, err] = simple_fixedpoint(phi0, p, N, maxit, tol, do_check_poincare);
			cpu_time = toc;
			data_fixp(i1, i2, 1) = data_fixp(i1, i2, 1) + err(end);
			data_fixp(i1, i2, 2) = data_fixp(i1, i2, 2) + length(err);
			data_fixp(i1, i2, 3) = data_fixp(i1, i2, 3) + cpu_time;
			fprintf("...done!\n")


			% end fixed point
			% Newton iteration for the coefficients of phi(z)
			fprintf("Computing coefficients of phi with Newton...")
			tic;
			[phi, err] = newton_iter(phi0, p, N, maxit, tol, do_check_poincare); 
			cpu_time = toc;
			data_newt(i1, i2, 1) = data_newt(i1, i2, 1) + err(end);
			data_newt(i1, i2, 2) = data_newt(i1, i2, 2) + length(err);
			data_newt(i1, i2, 3) = data_newt(i1, i2, 3) + cpu_time;
			fprintf("...done!\n")
			% end Newton
			
			fprintf("Computing coefficients of phi with the recursive formula...")
			tic;
			[phi, err] = recursive_solution_fft(p, N);
			cpu_time = toc;
			data_recu(i1, i2, 1) = data_recu(i1, i2, 1) + err;
			data_recu(i1, i2, 2) = data_recu(i1, i2, 2) + cpu_time;
			fprintf("...done!\n")

		end
		fprintf('---------Results for d = %d, and m = %d---------\n', d, m)
		fprintf('Fixed-point\n')
		fprintf('Avg residual = %.4e, \t Avg iterations = %.1f,\t Avg cpu time = %.3f s\n',...
		data_fixp(i1, i2, 1)/ntests, data_fixp(i1, i2, 2)/ntests, data_fixp(i1, i2, 3)/ntests);
		fprintf('Newton\n')
		fprintf('Avg residual = %.4e, \t Avg iterations = %.1f,\t Avg cpu time = %.3f s\n',...
		data_newt(i1, i2, 1)/ntests, data_newt(i1, i2, 2)/ntests, data_newt(i1, i2, 3)/ntests);
		fprintf('Recursive formula\n')
		fprintf('Avg residual = %.4e, \t Avg cpu time = %.3f s\n',...
		data_recu(i1, i2, 1)/ntests, data_recu(i1, i2, 2)/ntests);
		fprintf('------------------------------\n')
		i2 = i2 + 1;
	end
	i1 = i1 + 1;
end

data_fixp = data_fixp / ntests;
data_newt = data_newt / ntests;
data_recu = data_recu / ntests;


% Write fixed-point data in file
fileID = fopen('data/test_fix_vs_newton_timefx.txt','w');
fprintf(fileID, '& $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\ \n', avg_offsp);
for j = 1:size(data_fixp, 1)
	fprintf(fileID, '$%d$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\\\ \n', [degree(j), data_fixp(j, :, 3)]);
end
fclose(fileID);
fileID = fopen('data/test_fix_vs_newton_resfx.txt','w');
fprintf(fileID, '& $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\ \n', avg_offsp);
for j = 1:size(data_fixp, 1)
	fprintf(fileID, '$%d$ & $%1.0e$ & $%1.0e$ & $%1.0e$ & $%1.0e$ & $%1.0e$ \\\\ \n', [degree(j), data_fixp(j, :, 1)]);
end
fclose(fileID);
% Write Newton data in file
fileID = fopen('data/test_fix_vs_newton_timenw.txt','w');
fprintf(fileID, '& $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\ \n', avg_offsp);
for j = 1:size(data_fixp, 1)
	fprintf(fileID, '$%d$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\\\ \n', [degree(j), data_newt(j, :, 3)]);
end
fclose(fileID);
fileID = fopen('data/test_fix_vs_newton_resnw.txt','w');
fprintf(fileID, '& $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\ \n', avg_offsp);
for j = 1:size(data_fixp, 1)
	fprintf(fileID, '$%d$ & $%1.0e$ & $%1.0e$ & $%1.0e$ & $%1.0e$ & $%1.0e$ \\\\ \n', [degree(j), data_newt(j, :, 1)]);
end
fclose(fileID);
% Write recursive formula data in file
fileID = fopen('data/test_fix_vs_newton_timere.txt','w');
fprintf(fileID, '& $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\ \n', avg_offsp);
for j = 1:size(data_fixp, 1)
	fprintf(fileID, '$%d$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\\\ \n', [degree(j), data_newt(j, :, 2)]);
end
fclose(fileID);
fileID = fopen('data/test_fix_vs_newton_resre.txt','w');
fprintf(fileID, '& $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\ \n', avg_offsp);
for j = 1:size(data_fixp, 1)
	fprintf(fileID, '$%d$ & $%1.0e$ & $%1.0e$ & $%1.0e$ & $%1.0e$ & $%1.0e$ \\\\ \n', [degree(j), data_recu(j, :, 1)]);
end
fclose(fileID);



