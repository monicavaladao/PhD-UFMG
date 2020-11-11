function [] = launch_DEbest(problem, alg, seed, dir_output, filename)

% Get problem data
n = problem.n;
lb = problem.lb;
ub = problem.ub;
fobj = problem.fobj;
npop = problem.npop;
neval = problem.neval;
name = problem.name;

% Create initial sample
rng(seed, 'twister');
ssize = 4 * npop;
X = lhsdesign(ssize, n);
y = evaluation(fobj, X',lb,ub);

% Select initial population
[~,I] = sort(y,'ascend');
P = X(I(1:npop),:)';
jP = y(I(1:npop));

% Solve the problem
[XMIN,OBJMIN,G,GXMIN,GOBJMIN,GNfobj,numfobj,RUNTime] = algoritmo_de_best_normalizado(fobj,P,jP,lb(:),ub(:),npop,ssize,neval);


% Save results in a CSV file
if ~exist(dir_output, 'dir')
    mkdir(dir_output);
end

fid = fopen(strcat(dir_output, '/', filename), 'w+');
fprintf(fid, 'METAMODEL,PROB,NVAR,REP,NEVAL,ITER,BEST.OBJ,TOTAL.TIME.S\n');

history = struct();
history.iteractions = G;
history.best_x = GXMIN;
history.best_y = GOBJMIN;
history.neval = GNfobj;
history.debest_runtime = RUNTime;

for i = 1:length(GOBJMIN)
    fprintf(fid, '"%s","%s",%d,%d,%d,%d,%.6f,%.6f\n', ...
        alg, problem.name, n, seed, GNfobj(i), ...
        G(i), GOBJMIN(i), RUNTime(i));
end
fclose(fid);

% Save history in a MAT file
filename_mat = strrep(filename, '.csv', '.mat');
save(strcat(dir_output, '/', filename_mat), 'history');


end

