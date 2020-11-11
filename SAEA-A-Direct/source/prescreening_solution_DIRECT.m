function [x_b, y_b, pred_b, v_b, cr_b, imp_b, info] = prescreening_solution_DIRECT(fobj, P, P_V, CR, x_best, lb, ub, rule_choose_sample, sample_size, params, pool, info)
% PRESCREENING_SOLUTIONS: Select solution to evaluate on original function.
% Input:
%   fobj: Handle to the objective function
%   P : Offsprings
%	P_V: Behavior correction of each offspring in P
%	CR: Vector crossover rate (Nx1)
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
%	v_b: Behavior correction of x_b
%	cr_b: Crossover rate of x_b
%	imp_b: Percentage of improvment

% Identify parameters
L = params.LDS.l;                   % Dimension of the reduced low-dimensional space
type = params.metamodel;            % Type of metamodel
std_tol = params.tol_std;           % Threshold value of the variance

% Identfy the database
X = pool.X;
y = pool.y;

% To prevent a poorly metamodel sample
if info.pool.poorly_quality == 1
    sample_X = info.pool.Xsample;
    sample_y = info.pool.Ysample;
else    
    % Select the metamodel sample
    [sample_X, sample_y, info] = select_metamodel_sample(x_best, sample_size, std_tol, rule_choose_sample, params, info);
end


% Buil a Simplified Kriging metamodel
lb_d1 = -1; % To buil on one-dimension - the correct value is defined in build_metamodel_DIRECT
ub_d1 = 1;  % To buil on one-dimension - the correct value is defined in build_metamodel_DIRECT
[model_info] = build_metamodel_DIRECT(sample_X, sample_y, lb_d1, ub_d1, type, params);

% 
switch type
    case 'RBF_SRGTSToolbox'
        % Identfy the prediction function
        f_pred = model_info.fobjPredicao;

        % Evaluate P using prediction function
        [P_pred] = feval_all(f_pred, P);

        % Find solution with the lowest value of predicted function
        [min_sort, idx_sort] = sort(P_pred,'ascend'); 
        r_idx = randperm(3,1);
        idx_pred = idx_sort(r_idx);
        x_b = P(idx_pred,:);
        pred_b = P_pred(idx_pred);

        % Evaluate solution on the original function
        y_b = feval_all(fobj, x_b);
        
        % Behavior correction
        v_b = P_V(idx_pred,:);
        
        % Crossover rate
        cr_b = CR(idx_pred);

    otherwise % OK, UK1, UK2 and BK from ooDACE toolbox only.   
        % Identify the Expcted Improvement function
        f_EI = model_info.fobjEI;

        % Evaluate P using Expected Improvement function
        [P_EI, P_pred] = feval_all_two_output(f_EI, P);

        % Find solution with the highest value of Expected Improvement
        [max_EI, idx_EIsort] = sort(P_EI,'descend'); 
        r_idx = randperm(3,1);
        idx_EI = idx_EIsort(r_idx);
        x_b = P(idx_EI,:);
        pred_b = P_pred(idx_EI);

        % Evaluate solution on the original function
        y_b = feval_all(fobj, x_b);
        
        % Behavior correction
        v_b = P_V(idx_EI,:);
        
        % Crossover rate
        cr_b = CR(idx_EI);
end

%
fbest_initial = params.S_A.fbest_initial;
y_ref = min(fbest_initial,y_b);
imp_b = abs(((fbest_initial - y_ref)./fbest_initial).*100);
%imp_b = abs(((fbest_initial - y_b)./fbest_initial).*100);


end
