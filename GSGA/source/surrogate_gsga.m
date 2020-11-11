function [best_x, best_y, info] = surrogate_gsga(fobj, X, y, lb, ub, max_eval, varargin)
% SURROGATE_GSGA: Surrogate Assisted Evolutionary Algorithm  build over
% ooDACE toolbox and SRGTS toolbox.
%
% Input:
%   fobj: handle to the objective function
%   X: Sample (rows are entries and coluns are the variables)
%   y: Evaluation of each row in X
%   lb: Lower bounds
%   ub: Upper bounds
%   max_eval: Budget of objective function evaluations
%
%
% Output:
%   best_x: Best solution found.
%   best_y: Objective value of the best solution found.
%   info: A structure with additional information.
%
% References:
% [1] Liang2020

% Start timer
t0_start = cputime;

% Initialize structures used by the SAEA
[problem, params] = build_params_structure(fobj, X, y, lb, ub, max_eval, varargin{:});
info = build_info_structure(X, y, problem, params);

% Get some parameters
n = problem.n;          % Number of variables
N = params.pop_size;    % Population size of the EA
Ns = params.gsga.Ns;    % GSGA parameter
N1 = params.gsga.N1;    % GSGA parameter
N2 = params.gsga.N2;    % GSGA parameter
N3 = params.gsga.N3;    % GSGA parameter
pcross = params.gsga.pcross; % GSGA parameter
pmut = params.gsga.pmut; % GSGA parameter

% Select the population
[pop_X, pop_y] = select_population(X, y, N);
idx_best = 1;

% Find the current best solution
[value, idx] = min(y);
info.best_x = X(idx,:);
info.best_y = value;

% Update stats
info.history.iterations = [0];
info.history.best_x = [info.best_x];
info.history.best_y = [info.best_y];
info.history.neval = [info.neval];
info.history.mean_diff = [0];
info.history.metamodel_runtime = [0];
info.history.saea_runtime = [0];


% Logging
if params.verbose
    fprintf('-------------------------------------------------------------------------- \n');
    fprintf(' Iterations |      Best Obj. | Fun.Eval. |     Mean Diff. |    Runtime (s) \n');
    fprintf('-------------------------------------------------------------------------- \n')
    fprintf('% 11d | % 14.5f | % 9d |            --- | % 14.5f \n', 0, ...
        info.best_y, info.neval, (cputime - t0_start));
end

% Initialize counters
eval_counter = info.neval;
iter_counter = 1;

while eval_counter < params.max_eval     
    % Find the best solution in pop_X   
    [~,idx_sort] = sort(pop_y, 'ascend');
    pop_X = pop_X(idx_sort,:);
    pop_y = pop_y(idx_sort);
    best_x = pop_X(idx_sort(1),:);
    
    % Perform the Surrogate-Guided GA Updating Mechanism
    [RBF_X,RBF_y] = sga(pop_X,best_x,lb,ub,params,info);
    
    % Create the child population
    [P_child] = create_offsprings_ga(pop_X,RBF_X,lb,ub,n,N,N1,N2,pcross,pmut);
  
    % Perform the prescreening on Simplified Kriging (essa etapa retorna as novas soluções avaliadas na função original)
    [New_P,New_y,New_pred] = prescreening_simplified_kriging(fobj, lb, ub, P_child, params, info);
    
    % Update the evaluation counter
    info.neval = info.neval + Ns;
    eval_counter = eval_counter + Ns;
    
    % Update the pool of solutions
    info.pool = update_pool(pop_X, New_P, New_y, info.pool, params);

    % Update best solution
    aux_X = [New_P; info.best_x];
    aux_y = [New_y; info.best_y];
    [value, idx] = min(aux_y);
    info.best_x = aux_X(idx,:);
    info.best_y = value;
    
    % Update the history
    info.history.iterations = [info.history.iterations, iter_counter];
    info.history.neval = [info.history.neval, info.neval];
    info.history.best_x = [info.history.best_x; info.best_x];
    info.history.best_y = [info.history.best_y, info.best_y];
    info.history.mean_diff = [info.history.mean_diff, abs(mean(New_y(:) - New_pred(:)))];
    info.history.saea_runtime = [info.history.saea_runtime, cputime - t0_start];
    
    % Logging
    if params.verbose
        if info.history.best_y(end) < info.history.best_y(end-1)
            flag = '*';
        else
            flag = '';
        end
        fprintf('%1s% 10d | % 14.5f | % 9d | % 14.5f | % 14.5f \n', ...
            flag, iter_counter, info.best_y, info.neval, ...
            abs(mean(New_y(:) - New_pred(:))), (cputime - t0_start));
    end
        
    % Select the population of next generation
    [pop_X,pop_y] = select_gsga_population(pop_X,pop_y,New_P,New_y,N3);
    
    % Update iteration counter
    iter_counter = iter_counter + 1;
   

end

best_x = info.best_x;
best_y = info.best_y;

end