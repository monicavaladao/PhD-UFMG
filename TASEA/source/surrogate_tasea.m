function [best_x, best_y, info] = surrogate_tasea(fobj, X, y, lb, ub, max_eval, varargin)
% SURROGATE_TASEA: Surrogate Assisted Evolutionary Algorithm  build over
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
% Optional input (key/value pairs):
%   - EvolutionControl: Strategy used to control solutions.
%       Values: 'metamodel', 'random'.
%   - Metamodel: Type of metamodel.
%       Values: 'OrdinaryKriging', 'UniversalKriging1',
%       'UniversalKriging2', 'RBF'.
%   - Optimizer: Algorithm used to optimize the metamodel parameters. It is
%       used with Kriging metamodels only.
%       Values: 'sqp', 'fmincon', 'ga'.
%   - RBF: Type of RBF function. It is used with RBF metamodel only.
%       Values: 'Gaussian', 'GaussianCrossValidation', 'Multiquadric'.
%
% Output:
%   best_x: Best solution found.
%   best_y: Objective value of the best solution found.
%   info: A structure with additional information.
%
% References:
% [1] Yang, Z., Qiu, H., Gao L., Jiang, C. Zhang J.: Two-layer adaptive
%     surrogate-assisted evolutionary algorithm for high-dimensional
%     computationaly expensive problems. Journal of Global Optimization.
%     (74), 327 - 359(2019).


% Start timer
t0_start = cputime;

% Initialize structures used by the SAEA
[problem, params] = build_params_structure(fobj, X, y, lb, ub, max_eval, varargin{:});
info = build_info_structure(X, y, problem, params);

% Get some parameters
n = problem.n;          % Number of variables
N = params.pop_size;    % Population size of the EA
rule_choose_sample = params.choose_sample.rule; % Way to choose de metamodel sample

% Initial sample size
sample_size = params.sample_size;

% Sample size counter of local metamodel 
size_counter_ls = sample_size;

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
    best_x = pop_X(idx_best, :);
    %best_x = info.best_x;
    
    % Create offspring solutions using DE operators
    P = create_offsprings(pop_X, best_x, lb, ub, params.offsprings_per_solution, params.FM.Nf, params.FM.Nc, params.FM.Ngood);
    
    % Prescreening solution
    [x_b, y_b, pred_b, info] = prescreening_solution(fobj, P, best_x, lb, ub, rule_choose_sample, sample_size, params, info.pool, info);
    
    % Update the evaluation counter
    info.neval = info.neval + 1;
    eval_counter = eval_counter + 1;
       
    % Update the pool of solutions
    info.pool = update_pool(pop_X, x_b, y_b, info.pool, params);

    % Update best solution
    aux_X = [x_b; info.best_x];
    aux_y = [y_b; info.best_y];
    [value, idx] = min(aux_y);
    info.best_x = aux_X(idx,:);
    info.best_y = value;
    
    % Update the Feedback Mechanism parameters
    if idx == 1
        params.FM.Nf = 1;
        params.FM.Ngood = 1;
        rule_choose_sample = params.choose_sample.rule_local_metamodel;
        size_counter_ls = size_counter_ls + 1;  
        sample_size = min(size_counter_ls,params.sample_size_local_metamodel);
    else
        params.FM.Nf = params.FM.Nf + 1;
        params.FM.Ngood = 0;
        rule_choose_sample = params.choose_sample.rule;
        sample_size = params.sample_size;
    end

    % Update the history
    info.history.iterations = [info.history.iterations, iter_counter];
    info.history.neval = [info.history.neval, info.neval];
    info.history.best_x = [info.history.best_x; info.best_x];
    info.history.best_y = [info.history.best_y, info.best_y];
    info.history.mean_diff = [info.history.mean_diff, abs(pred_b - y_b)];
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
            abs(pred_b - y_b), (cputime - t0_start));
    end
    
    % Select the population
    [pop_X, pop_y] = select_population(info.pool.X, info.pool.y, N);
    idx_best = 1;
    
    % Update iteration counter
    iter_counter = iter_counter + 1;

     
end

best_x = info.best_x;
best_y = info.best_y;

end