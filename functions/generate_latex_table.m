function generate_latex_table(data, file_path)
	if size(data, 2) < 26
		data = [data, nan(size(data, 1), 26 - size(data, 2))];
	end
	fileID = fopen(file_path,'w');
	for j = 1:size(data, 1)
		fprintf(fileID, '& $%d$ & $%.2f$ & $%.1e$ & $%d$ & $%.2f$ & $%.1e$ & $%.2f$ & $%.1e$ & $%d$ & $%.1f$ & $%d$ & $%.2f$ & $%.1e$ & $%d$ & $%.1f$ & $%d$ & $%.2f$ & $%.1e$ & $%d$ & $%.1f$ & $%d$ & $%.2f$ & $%.1e$ & $%d$ & $%.1f$ & $%d$ \\\\ \n', data(j, :));
	end
	fclose(fileID);
end
