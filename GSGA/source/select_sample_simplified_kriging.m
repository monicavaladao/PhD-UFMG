function[sample_X, sample_y, info] = select_sample_simplified_kriging(sample_size, std_tol, params, info)
% SELECT_SAMPLE_SIMPLIFIED_KRIGING: Select the metamodel sample.
% Input: 
%   sample_size: Metamodel sample size.
%   std_tol: Threshold value of the variance
%   info: A structure used to keep the progress of the GSGA algorithm
% Output:
%   sample_X: The sample selected to build the metamodel (rows are entries 
%       and coluns are the variables).
%   sample_y: Evaluate of each row in X.
%   info: A structure used to keep the progress of the GSGA algorithm

% Identify the pool of solutions
pool = info.pool;

% Select a sample to build/update the metamodel
[sample_X, sample_y] = choose_newest(pool, sample_size, std_tol);


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
% Auxiliar function
% -------------------------------------------------------------------------
%
%
function [X,y] = choose_newest(pool,sample_size, std_tol)
% CHOOSE_NEWESTS: Choose a sample to create/update the metamodel.
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
