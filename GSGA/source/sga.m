function [Pp,FPp] = sga(P,x_best,lb,ub,params,info)
% SGA: Perform the Surrogate-Guided GA Updating Mechanism.
%Input:
%   P: The current population
%   x_best: The best soltion in P
%   lb: Lower bound
%   ub: Upper bound
%   info: A structure with informations
% Output:
%   Pp: The predicted population
%   FPp: The predicted value of each solution in Pp

% Identify the RBF metamodel
type = params.metamodel_rbf;   

% Size of population
[N,n] = size(P);

% Size of the metamodel sample
%sample_size = params.sample_size.rbf;
pool_size = length(info.pool.y);
sample_size = min(5*n + 10,pool_size);

% Ver melhor a definição dos parâmtros do fmincon ou se utiliza outro
% algoritmo para essa otimização.
%options = optimoptions('fmincon','Display','off', 'maxIterations', 100,...
%    'maxFunEvals', 100, 'functionTolerance', 1e-4,...
%    'gradobj', 'off', 'algorithm','active-set','diagnostics', 'off','derivativecheck', 'off');

for i = 1:N
    % Identify the solution x
    x = P(i,:);
    
    % Identify the farthest solution from x
    d_dist = sqrt(sum((repmat(x, N, 1) - P) .^ 2, 2));
    [d_sort, index_sort] = sort(d_dist,'descend');
    D_max = d_sort(1);
    r =  0.5*(D_max/(sqrt(n)*nthroot(N - 1, n)));
    lb_xr = x_best - r*ones(1,n);
    ub_xr = x_best + r*ones(1,n);
%    lb_xr = x - r*ones(1,n);
%    ub_xr = x + r*ones(1,n);

    % Truncation
    lb_xr(lb_xr < lb) = lb(lb_xr < lb);
    ub_xr(ub_xr > ub) = ub(ub_xr > ub);
    
    % Select metamodel sample
    [sample_X, sample_y] = select_sample_rbf(x, lb_xr, ub_xr, info, params, sample_size);
    
    % Buil a RBF metamodel
    [model_info] = build_metamodel_DIRECT(sample_X, sample_y, lb_xr, ub_xr, type, params);
    
    % Identfy the prediction function
    f_pred = model_info.fobjPredicao;
    
    % Optimize the prediction function
%    [x_m,f_m] = fmincon(f_pred,x_best,[],[],[],[],lb_xr,ub_xr,[],options);
    
    [x_m,f_m] = sl_pso_optimize(f_pred,lb_xr,ub_xr);
    
    % Store the solutions
    Pp(i,:) = x_m(:)';
    FPp(i,1) = f_m;
    
end

end