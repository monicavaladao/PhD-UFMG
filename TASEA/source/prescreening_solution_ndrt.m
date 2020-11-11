% % function [x_b, y_b, pred_b] = prescreening_solution(fobj, P, x_best, lb, ub, rule_choose_sample, sample_size, params, pool, info)
% % % PRESCREENING_SOLUTIONS: Select solution to evaluate on original function.
% % % Input:
% % %   fobj: Handle to the objective function
% % %   P : Offsprings
% % %   x_best: Best solution in the current population
% % %   lb: Lower bounds
% % %   ub: Upper bounds
% % %   N: Number of solutions in the current population
% % %   n: Number of variables
% % %   rule_choose_sample: Rule to choose sample
% % %   sample_size: Metamodel sample size
% % %   params: A structure with fields that keep parameters used by the SAEA
% % %       (including toolbox ooDACE)
% % %   info: A structure used to keep the progress of the SAEA algorithm
% % % Output:
% % %   x_b: Solution whith best predicted value
% % %   y_b: Evaluation of x_b on objective function fobj
% % %   pred_b: Evaluation of x_b on metamodel
% % 
% % % Identify parameters
% % L = params.LDS.l;                   % Dimension of the reduced low-dimensional space
% % type = params.metamodel;            % Type of metamodel
% % std_tol = params.tol_std;           % Threshold value of the variance
% % 
% % % Identfy the database
% % X = pool.X;
% % y = pool.y;
% % 
% % % Size of database
% % [Nx,~] = size(X);
% % 
% % % Apply the dimensional reduction technique
% % info_dr = struct();
% % info_dr = info;
% % Aux_XP = [X;P;x_best;lb;ub];
% % %S = sammon(Aux_XP, L, 10, 'seconds')'; 
% % %S = sammon(Aux_XP, L, 5, 'seconds')'; 
% % S = sammon(Aux_XP, L, 100, 'steps')';
% % %S = sammon(Aux_XP, L, 1e-3,'errlimit')'; 
% % %S = sammon(Aux_XP, L, 1e-5,'errchange')'; 
% % 
% % 
% % S = S';
% % X_dr = S(1:Nx,:);
% % P_dr = S((Nx+1):(end-3),:);
% % x_best_dr = S((end - 2),:);
% % lb_dr = S((end-1),:);
% % ub_dr = S(end,:);
% % info_dr.pool.X = X_dr;
% % 
% % % Select the metamodel sample
% % [sample_X, sample_y] = select_metamodel_sample(x_best_dr, sample_size, std_tol, rule_choose_sample, info_dr);
% % 
% % % Build a metamodel
% % [model_info] = build_metamodel(sample_X, sample_y, lb_dr, ub_dr, type, params);
% % 
% % % Identify the Expcted Improvement function
% % f_EI = model_info.fobjEI;
% %  
% % % Evaluate P_dr using Expected Improvement function
% % offs_X = P_dr;
% % [offs_EI, offs_pred] = feval_all_two_output(f_EI, offs_X);
% % 
% % % Find solution with the highest value of Expected Improvement
% % [max_EI, idx_EIsort] = sort(offs_EI,'descend'); 
% % r_idx = randperm(3,1);
% % idx_EI = idx_EIsort(r_idx);
% % x_dr = offs_X(idx_EI,:);
% % 
% % % Corresponding x_dr in the original space
% % x_b = P(idx_EI,:);
% % pred_b = offs_pred(idx_EI);
% % 
% % % Evaluate solution on the original function
% % y_b = feval_all(fobj, x_b);
% % 
% % end

function [x_b, y_b, pred_b] = prescreening_solution_ndrt(fobj, P, x_best, lb, ub, rule_choose_sample, sample_size, params, pool, info)
% PRESCREENING_SOLUTIONS: Select solution to evaluate on original function.
% Input:
%   fobj: Handle to the objective function
%   P : Offsprings
%   x_best: Best solution in the current population
%   lb: Lower bounds
%   ub: Upper bounds
%   N: Number of solutions in the current population
%   n: Number of variables
%   rule_choose_sample: Rule to choose sample
%   sample_size: Metamodel sample size
%   params: A structure with fields that keep parameters used by the SAEA
%       (including toolbox ooDACE)
%   info: A structure used to keep the progress of the SAEA algorithm
% Output:
%   x_b: Solution whith best predicted value
%   y_b: Evaluation of x_b on objective function fobj
%   pred_b: Evaluation of x_b on metamodel

% Identify parameters
L = params.LDS.l;                   % Dimension of the reduced low-dimensional space
type = params.metamodel;            % Type of metamodel
std_tol = params.tol_std;           % Threshold value of the variance

% Identfy the database
X = pool.X;
y = pool.y;

% Size of database
%[Nx,~] = size(X);

% Select the metamodel sample
[sample_X, sample_y] = select_metamodel_sample(x_best, sample_size, std_tol, rule_choose_sample, info);

% Build a metamodel
[model_info] = build_metamodel(sample_X, sample_y, lb, ub, type, params);

% 
switch type
    case 'RBF_SRGTSToolbox'
        % Identfy the prediction function
        f_pred = model_info.fobjPredicao;

        % Evaluate P using prediction function
        offs_X = P;
        [offs_pred] = feval_all(f_pred, offs_X);

        % Find solution with the lowest value of predicted function
        [min_sort, idx_sort] = sort(offs_pred,'ascend'); 
        r_idx = randperm(3,1);
        idx_pred = idx_sort(r_idx);
        x_b = offs_X(idx_pred,:);
        pred_b = offs_pred(idx_pred);

        % Evaluate solution on the original function
        y_b = feval_all(fobj, x_b);

    otherwise % OK, UK1, UK2 and BK from ooDACE toolbxonly.   
        % Identify the Expcted Improvement function
        f_EI = model_info.fobjEI;

        % Evaluate P using Expected Improvement function
        offs_X = P;
        [offs_EI, offs_pred] = feval_all_two_output(f_EI, offs_X);

        % Find solution with the highest value of Expected Improvement
        [max_EI, idx_EIsort] = sort(offs_EI,'descend'); 
        r_idx = randperm(3,1);
        idx_EI = idx_EIsort(r_idx);
        x_b = offs_X(idx_EI,:);
        pred_b = offs_pred(idx_EI);

        % Evaluate solution on the original function
        y_b = feval_all(fobj, x_b);
end

end
