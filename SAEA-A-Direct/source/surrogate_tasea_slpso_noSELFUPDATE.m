function [best_x, best_y, info] = surrogate_tasea_slpso_noSELFUPDATE(fobj, X, y, lb, ub, max_eval, varargin)
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
% Implement the SAEA_SL-PSO_DE-SELF_UPDATE


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
V = info.pool.behavior_vector;
CR = info.pool.CR;
[pop_X, pop_y, pop_V, pop_CR] = select_population(X, y, V, CR, N);
idx_best = 1;

% Find the current best solution
[value, idx] = min(y);
info.best_x = X(idx,:);
info.best_y = value;

% Store information used in self update parameters
params.CR = pop_CR;
params.S_A.CRrec = pop_CR;
params.S_A.fbest_initial = info.best_y;
params.S_A.Pimp = abs(((repmat(info.best_y,N,1) - pop_y(:))./info.best_y).*100);

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
eval_counter = info.neval; % Evaluation counter
iter_counter = 1;          % Iteration
DE_counter = 0;            % Used in DE operators
any_DE_used = 0;           % Used in DE operators
bestever = 1e200;          % Used in SL_PSO operator
while eval_counter < params.max_eval
       
    % Find the best solution in pop_X
    best_x = pop_X(idx_best, :);
    %best_x = info.best_x;
    
    % Create offspring solutions using SL_PSO and operators
    [P,P_V,P_CR,Nf,DE_counter, any_DE_used] = create_offsprings_noSELFUPDATE(pop_X, pop_y, pop_V, pop_CR, best_x, lb, ub, ...
                                 DE_counter, bestever, params.offsprings_per_solution, params.FM.Nf, params.FM.Nc, params.FM.Ngood, any_DE_used);
    params.FM.Nf = Nf;
    params.CR = P_CR;
        
    [x_b, y_b, pred_b, v_b, cr_b, imp_b, info] = prescreening_solution_DIRECT(fobj, P, P_V, P_CR, best_x,...
                                           lb, ub, rule_choose_sample, sample_size, params, info.pool, info);
    
    % Update the evaluation counter
    info.neval = info.neval + 1;
    eval_counter = eval_counter + 1;
       

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
        
%         %%%%
%         params.S_A.CRrec = [params.S_A.CRrec;cr_b];
%         params.S_A.Pimp = [params.S_A.Pimp;imp_b];
    else
        params.FM.Nf = params.FM.Nf + 1;
        params.FM.Ngood = 0;
        rule_choose_sample = params.choose_sample.rule;
        sample_size = params.sample_size;
    end
    
    % Self update parameters
% %     if any_DE_used == 1
% %         params.S_A.CRrec = [params.S_A.CRrec;cr_b];
% %         params.S_A.fbest_initial = info.best_y;
% %         params.S_A.Pimp = [params.S_A.Pimp;imp_b];
% %         [CR, CRm, S_A] = crossover_rate_self_adaptation(params.CR, params.CRm, params.S_A, N, n, DE_counter);
% %         params.CR = CR;
% %         params.CRm = CRm;
% %         params.S_A.CRrec = S_A.CRrec;
% %         params.S_A.Pimp = S_A.Pimp;
% %         
% %         any_DE_used = 0;
% %     else
% %         params.S_A.CRrec = [params.S_A.CRrec;cr_b];
% %         params.S_A.fbest_initial = info.best_y;
% %         params.S_A.Pimp = [params.S_A.Pimp;imp_b];
% %         CR = params.CR;
% %     end
        
    
    %%%
     params.S_A.CRrec = [params.S_A.CRrec;cr_b];
     params.S_A.fbest_initial = info.best_y;
     params.S_A.Pimp = [params.S_A.Pimp;imp_b];
%    [CR, CRm, S_A] = crossover_rate_self_adaptation(params.CR, params.CRm, params.S_A, N, n, DE_counter);
    [CR, CRm, S_A] = crossover_rate_self_adaptation(params.CR, params.CRm, params.S_A, N, n, iter_counter);
    params.CR = CR;
    params.CRm = CRm;
    params.S_A.CRrec = S_A.CRrec;
    params.S_A.Pimp = S_A.Pimp;    
    %%%
      
    
    % Update the pool of solutions
    pop_CR = CR;
    info.pool = update_pool(pop_X, pop_V, pop_CR, x_b, y_b, v_b, cr_b, info.pool, params);

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
    [pop_X, pop_y, pop_V, pop_CR] = select_population(info.pool.X, info.pool.y, info.pool.behavior_vector, info.pool.CR, N);
    idx_best = 1;
    
    bestever = info.best_y;
    
    % Update iteration counter
    iter_counter = iter_counter + 1;

     
end

best_x = info.best_x;
best_y = info.best_y;

end