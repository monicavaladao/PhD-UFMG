function [New_P,New_y,New_pred] = prescreening_simplified_kriging(fobj, lb, ub, P_child, params,info)
% PRESCREENING_SIMPLIFIED_KRIGING: Perform the Prescreening Strategy Based on Simplified Kriging
% Input: 
%   fobj: The true function
%   lb: Lower bound
%   ub: Upper bound
%   P_child: The child population
% Output:
%   New_P: The new solutions evalauted on fobj
%   New_y: The fobj value of each solution in New_P
%   New_pred: The predicted value of each solution in New_P

% Identify the RBF metamodel
type = params.metamodel_kriging;   

% Number of solutions to be evaluated on the original function
Ns = params.gsga.Ns;

% Threshold value of the variance
std_tol = params.tol_std;  

% Size of the metamodel sample
%sample_size = params.sample_size.kriging;
pool_size = length(info.pool.y);
sample_size = pool_size;

% Select the metamodel sample
[sample_X, sample_y, info] = select_sample_simplified_kriging(sample_size, std_tol, params, info);

% Buil a Simplified Kriging metamodel
lb_d1 = -1; % To buil on one-dimension - the correct value is defined in build_metamodel_DIRECT
ub_d1 = 1;  % To buil on one-dimension - the correct value is defined in build_metamodel_DIRECT
[model_info] = build_metamodel_DIRECT(sample_X, sample_y, lb_d1, ub_d1, type, params);

% Identify the Expcted Improvement function
f_EI = model_info.fobjEI;

% Evaluate P using Expected Improvement function
[child_EI, child_pred] = feval_all_two_output(f_EI, P_child);

% Find solution with the highest value of Expected Improvement
[max_EI, idx_EIsort] = sort(child_EI,'descend');

% Identify solutions to be evaluated on the original function
idx_Ns = idx_EIsort(1:Ns);
New_P = P_child(idx_Ns,:);
New_pred = child_pred(idx_Ns);

% Evauate the new solutions on the original function
New_y = feval_all(fobj, New_P);

end