function[sample_X, sample_y, info] = select_metamodel_sample(x_best, sample_size, std_tol, rule_choose_sample, params, info)
% SELECT_METAMODEL_SAMPLE: Select the metamodel sample.
% Input: 
%   x_best: Best solution in current populaton of EA
%   sample_size: Metamodel sample size.
%   std_tol: Threshold value of the variance
%   rule_choose_sample: Rule to choose sample
%   info: A structure used to keep the progress of the SAEA algorithm
% Output:
%   sample_X: The sample selected to build the metamodel (rows are entries 
%       and coluns are the variables).
%   sample_y: Evaluate of each row in X.

% Identify the pool of solutions
pool = info.pool;

% Select a sample to build/update the metamodel
switch rule_choose_sample
    case 'nearest'
        [sample_X, sample_y] = choose_nearest(x_best, pool, sample_size, std_tol);
    case 'newest'
        [sample_X, sample_y] = choose_newest(pool, sample_size, std_tol);
end

% To prevent a poorly metamodel sample
s_size = length(sample_y);
if s_size >= params.sample_max_sample_size
    info.pool.poorly_quality = 1;
else
    info.pool.Xsample = sample_X;
    info.pool.Ysample = sample_y;
end

end

% -------------------------------------------------------------------------
% Auxiliar functions
% -------------------------------------------------------------------------
%
%
function [X,y] = choose_newest(pool,sample_size, std_tol)
% CHOOSE_METAMODEL_SAMPLE: Choose a sample to create/update the metamodel.
% This function chooses the newest solutions into the pool to compose the
% sample.
%
% Input: 
%   pool: Pool of solutions.
%   sample_size: Metamodel sample size
%   std_tol: Threshold value of the variance
% 
% Output:
%   X: The sample selected to build the metamodel (rows are entries 
%       and coluns are the variables).
%   y: Evaluate of each row in X.
% 

% Sort solutions by their age
[~,idx] = sort(pool.age,'descend');
auxpool_X = pool.X(idx,:);
auxpool_y = pool.y(idx);
 
% Select the sample_size newest solution
X = auxpool_X(1:sample_size,:);
y = auxpool_y(1:sample_size);
 
% Number of solutions in the pool
pool_size = size(auxpool_X,1);

% If standart deviation of X or y is zero, then the metamodel sample must be redefined.
while (any(std(X) < std_tol) || any(isnan(std(X))) || std(y) < std_tol || isnan(std(y))) && sample_size < pool_size
    sample_size = sample_size + 1;
    X = auxpool_X(1:sample_size,:);
    y = auxpool_y(1:sample_size);
end

end
%
function [X,y] = choose_nearest(x_best, pool, sample_size, std_tol)
% CHOOSE_NEAREST: Choose a sample to create/update the metamodel.
% This function chooses the nearest solutions of x_best
%
% Input: 
%   pop_X: Current population of EA
%   pool: Pool of solutions
%   in the population
%   sample_size: Metamodel sample size
%   std_tol: Threshold value of the variance
% 
% Output:
%   X: The sample selected to build the metamodel (rows are entries 
%       and coluns are the variables).
%   y: Evaluate of each row in X.
% 

% Pool size N_pool
[N_pool, n] = size(pool.X);


% % Distance between x_best and pool
% d_xbest_pool = zeros(N_pool, 1);

% Identify the nearest solutions
d_xbest_pool = sqrt(sum((repmat(x_best, N_pool, 1) - pool.X) .^ 2, 2));
[~, index_sort] = sort(d_xbest_pool,'ascend');

% Correct the sample size
sample_size = min(sample_size,N_pool);

% Identify index
index_sample = index_sort(1:sample_size);

% Stores in X the solutions identified by index_sample
X = pool.X(index_sample,:);
y = pool.y(index_sample);

% Number of solutions in the pool
auxpool_X = pool.X;
auxpool_y = pool.y;
auxpool_X(index_sample,:) = [];
auxpool_y(index_sample) = [];
auxpool_X = [X;auxpool_X];
auxpool_y = [y;auxpool_y];
pool_size = size(auxpool_X,1);

sample_size = length(index_sample);
% If standart deviation of X or y is zero, then the metamodel sample must be redefined.
while (any(std(X) < std_tol) || any(isnan(std(X))) || std(y) < std_tol || isnan(std(y))) && sample_size < pool_size
    sample_size = sample_size + 1;
    X = auxpool_X(1:sample_size,:);
    y = auxpool_y(1:sample_size);
end


end
