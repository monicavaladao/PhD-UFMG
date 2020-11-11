% -------------------------------------------------------------------------
% Get all CSV files into a single one
fprintf('Joining CSV files into a single file...\n');

% Output directory in which results are saved
dir_output = './results_10_30/analytic';


% List of files with individual results
files = dir(strcat(dir_output, '/all/*.csv'));

% Write file with all results
fid = fopen(strcat(dir_output, '/results_10_30.csv'), 'w+');
fprintf(fid, 'METAMODEL,PROB,NVAR,REP,NEVAL,ITER,BEST.OBJ,MEAN.DIFF,TOTAL.TIME.S\n');

for i = 1:size(files, 1)
    filename = fullfile(files(i).folder, files(i).name);
    if exist(filename, 'file')
        cfid = fopen(filename, 'r');
        cline = fgetl(cfid); % ignore CSV header
        cline = fgetl(cfid); % first line after header
        while ischar(cline) && ~isempty(cline)
            fprintf(fid, cline); % write line
            fprintf(fid, '\n');
            cline = fgetl(cfid); % read next line
        end
        fclose(cfid);
    end
end

fclose(fid);
