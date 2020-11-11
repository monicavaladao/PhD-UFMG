function [X, y] = select_sample_rbf(x, lb_xr, ub_xr, info, params, sample_size)
% SELECT_SAMPLE_RBF: Select the sample of the RBF metamodel.

% Pool of solutions
auxpool_X = info.pool.X;
auxpool_y = info.pool.y;

% Number of solutions in the pool
pool_size = length(auxpool_y);

% Identify the nearest solutions
dist_x = sqrt(sum((repmat(x, pool_size, 1) - auxpool_X) .^ 2, 2));
[~, index_sort] = sort(dist_x,'ascend');
auxpool_X = auxpool_X(index_sort,:);
auxpool_y = auxpool_y(index_sort);

% Threshold value of the variance
std_tol = params.tol_std;  

% Auxiliary matrix
Aux_X = auxpool_X;
Aux_y = auxpool_y;

idx = [];
for i = 1:pool_size
    if any(auxpool_X(i,:) < lb_xr) || any(auxpool_X(i,:) > ub_xr)
        idx = [idx;i];
    end
end

Aux_X(idx,:) = [];
Aux_y(idx) = [];
X = Aux_X;
y = Aux_y;
N_y = length(y);

if N_y < sample_size
    X = auxpool_X(1:sample_size,:);
    y = auxpool_y(1:sample_size);
end

% If standart deviation of X or y is zero, then the metamodel sample must be redefined.
while (any(std(X) < std_tol) || any(isnan(std(X))) || std(y) < std_tol || isnan(std(y))) && sample_size < pool_size
    sample_size = sample_size + 1;
    X = auxpool_X(1:sample_size,:);
    y = auxpool_y(1:sample_size);
end


end