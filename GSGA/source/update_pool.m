function [pool] = update_pool(pop_X, chosen_X, chosen_y, pool, params)
% UPDADE_POOL: Upadate the pool of solutions.
%
% Input:
%   chosen_X: New solutions
%   chosen_y: Evaluate of each row in chosen_X
%   pool: Pool of solutions.


% Identify the rule to choose sample
rule_choose_sample = params.choose_sample.rule;

% Update the pool
switch rule_choose_sample
    case 'kmeans'
        [pool] = keeping_nearest_to_centroid(pop_X, chosen_X, chosen_y, pool, params);
    case 'lowest'
        [pool] = keeping_lowest(chosen_X, chosen_y, pool, params);
    case 'nearest'
        [pool] = keeping_nearest_to_pop(pop_X, chosen_X, chosen_y, pool, params); 
    case 'newest'
        [pool] = keeping_newest(chosen_X, chosen_y, pool, params);
end


end

% -------------------------------------------------------------------------
% Auxiliar functions
% -------------------------------------------------------------------------
%
function status = is_pool_ok(pool, tol_std)
    status = ~(any(std(pool.X) < tol_std) || any(isnan(std(pool.X))) || std(pool.y) < tol_std || isnan(std(pool.y)));
end
%
function [pool] = keeping_nearest_to_centroid(pop_X, chosen_X, chosen_y, pool, params)
% KEEPING_NEAREST_TO_CENTROID: Upadate the pool of solutions.
% This function updates the pool by keeping the solutions that are nearest 
% to centroid of the current EA population.

% Pool size N_pool and pop_X size N_pop
[N_pool, n] = size(pool.X);
[N_pop, ~] = size(pop_X);

% Remove repeated solutions from chosen_X
[chosen_X, idx, ~] = unique(chosen_X, 'rows', 'stable');
chosen_y = chosen_y(idx);
[N_chosen, ~] = size(chosen_X);

% Calculate the centroid
c = mean(pop_X,1);

% Distance between pool and c, for each solution in pool
D_pool = zeros(N_pool, 1);
for i = 1:N_pool
    D_pool(i) = sqrt(sum(pool.X(i,:) - c).^ 2);
end
[D_pool,idb_pool] = sort(D_pool,'ascend');
pool.X = pool.X(idb_pool,:);
pool.y = pool.y(idb_pool);

% Distance between chosen_X and c, for each solution in chosen_X
D_chosen = zeros(N_chosen, 1);
for i = 1:N_chosen
    D_chosen(i) = sqrt(sum(chosen_X(i,:) - c).^ 2);
end
[D_chosen,idb_chosen] = sort(D_chosen,'ascend');
chosen_X = chosen_X(idb_chosen,:);
chosen_y = chosen_y(idb_chosen);

% Try to insert solutions in chosen_X into pool
has_pool_changed = 0;
for i = 1:N_chosen
    
    % Distance between chosen_X(i,:) and all solutions in pool
    [d, idx] = min(sqrt(sum((repmat(chosen_X(i,:), N_pool, 1) - pool.X) .^ 2, 2)));
    
    if d <= params.tol_ratio
        
        % Check if chosen_y(i) <= pool.y(idx)
        if  chosen_y(i) <= pool.y(idx)
        
            % Keep a copy to undo the replacement, if necessary
            bkp_x = pool.X(idx,:);
            bkp_y = pool.y(idx);
            bkp_age = pool.age(idx);
            bkp_D_pool = D_pool(idx); 

            % Replace the solution in the pool
            pool.X(idx,:) = chosen_X(i,:);
            pool.y(idx) = chosen_y(i);
            pool.age(idx) = max(pool.age) + 1;
            D_pool(idx) = D_chosen(i);

            % Check the replacement
            if ~is_pool_ok(pool, params.tol_std)
                pool.X(idx,:) = bkp_x;
                pool.y(idx) = bkp_y;
                pool.age(idx) = bkp_age;
                D_pool(idx) = bkp_D_pool;
            else 
                has_pool_changed = 1;
            end
            
        end

    else
        
        if N_pool < params.max_pool_size
            
            % Add the solution into the pool
            idx = N_pool + 1;
            pool.X(idx, :) = chosen_X(i,:);
            pool.y(idx) = chosen_y(i);
            pool.age(idx) = max(pool.age) + 1;
            D_pool(idx) = D_chosen(i);
            
            
            % Check insertion
            if ~is_pool_ok(pool, params.tol_std)
                pool.X(idx,:) = [];
                pool.y(idx) = [];
                pool.age(idx) = [];
                D_pool(idx) = [];
            else
                N_pool = N_pool + 1;
                has_pool_changed = 1;
            end
            
        else
             % Find the solution in the pool that is the farthest from
             % centroid of pop_X
            [value, idx] = max(D_pool);
            
            % Keep a copy to undo the replacement, if necessary
            bkp_x = pool.X(idx,:);
            bkp_y = pool.y(idx);
            bkp_age = pool.age(idx);
            bkp_D_pool = value;

            % Replace the solution in the pool
            pool.X(idx,:) = chosen_X(i,:);
            pool.y(idx) = chosen_y(i);
            pool.age(idx) = max(pool.age) + 1;
            D_pool(idx) = D_chosen(i);
            
            % Check the replacement
            if ~is_pool_ok(pool, params.tol_std)
                pool.X(idx,:) = bkp_x;
                pool.y(idx) = bkp_y;
                pool.age(idx) = bkp_age;
                D_pool(idx) = bkp_D_pool;
            else 
                has_pool_changed = 1;
            end
         
        end
        
    end
end

% Update age status in pool structure
if has_pool_changed
    pool.next_age = max(pool.age) + 1;
end

end
%
function [pool] = keeping_lowest(chosen_X, chosen_y, pool, params)
% KEEPING_NEWEST: Upadate the pool of solutions.
% This function updates the pool by keeping the best solutions (i.e, the
% solutions with the lowest values).

% Pool size N_pool
[N_pool, n] = size(pool.X);

% Remove repeated solutions from chosen_X
[chosen_X, idx, ~] = unique(chosen_X, 'rows', 'stable');
chosen_y = chosen_y(idx);
[N_chosen, ~] = size(chosen_X);

% Try to insert solutions in chosen_X into pool
has_pool_changed = 0;
for i = 1:N_chosen
    
    % Distance between chosen_X(i,:) and all solutions in pool
    [d, idx] = min(sqrt(sum((repmat(chosen_X(i,:), N_pool, 1) - pool.X) .^ 2, 2)));
    
    if d <= params.tol_ratio
        
        % Check if chosen_y(i) <= pool.y(idx)
        if  chosen_y(i) <= pool.y(idx)
        
            % Keep a copy to undo the replacement, if necessary
            bkp_x = pool.X(idx,:);
            bkp_y = pool.y(idx);
            bkp_age = pool.age(idx);

            % Replace the solution in the pool
            pool.X(idx,:) = chosen_X(i,:);
            pool.y(idx) = chosen_y(i);
            pool.age(idx) = max(pool.age) + 1;

            % Check the replacement
            if ~is_pool_ok(pool, params.tol_std)
                pool.X(idx,:) = bkp_x;
                pool.y(idx) = bkp_y;
                pool.age(idx) = bkp_age;
            else 
                has_pool_changed = 1;
            end
            
        end

    else
        
        if N_pool < params.max_pool_size
            
            % Add the solution into the pool
            idx = N_pool + 1;
            pool.X(idx, :) = chosen_X(i,:);
            pool.y(idx) = chosen_y(i);
            pool.age(idx) = max(pool.age) + 1;
            
            % Check insertion
            if ~is_pool_ok(pool, params.tol_std)
                pool.X(idx,:) = [];
                pool.y(idx) = [];
                pool.age(idx) = [];
            else
                N_pool = N_pool + 1;
                has_pool_changed = 1;
            end
            
        else
            % Find the worst solution in pool
            [value, idx] = max(pool.y);
            
            % Keep a copy to undo the replacement, if necessary
            bkp_x = pool.X(idx,:);
            bkp_y = pool.y(idx);
            bkp_age = pool.age(idx);

            % Replace the solution in the pool
            pool.X(idx,:) = chosen_X(i,:);
            pool.y(idx) = chosen_y(i);
            pool.age(idx) = max(pool.age) + 1;

            % Check the replacement
            if ~is_pool_ok(pool, params.tol_std)
                pool.X(idx,:) = bkp_x;
                pool.y(idx) = bkp_y;
                pool.age(idx) = bkp_age;
            else 
                has_pool_changed = 1;
            end
         
        end
        
    end
end

% Update age status in pool structure
if has_pool_changed
    pool.next_age = max(pool.age) + 1;
end

end
%
function [pool] = keeping_nearest_to_pop(pop_X, chosen_X, chosen_y, pool, params)
% KEEPING_NEAREST_TO_POP: Upadate the pool of solutions.
% This function updates the pool by keeping the solutions that are nearest
% to the current population of EA.

% Pool size N_pool and pop_X size N_pop
[N_pool, n] = size(pool.X);
[N_pop, ~] = size(pop_X);

% Distance between pool and pop_X, for each solution in pool
d_pool_pop = zeros(N_pool, 1);
for i = 1:N_pool
    d_pool_pop(i) = min(sqrt(sum((repmat(pool.X(i,:), N_pop, 1) - pop_X) .^ 2, 2)));
end

% Remove repeated solutions from chosen_X
[chosen_X, idx, ~] = unique(chosen_X, 'rows', 'stable');
chosen_y = chosen_y(idx);
[N_chosen, ~] = size(chosen_X);

% Try to insert solutions in chosen_X into pool
has_pool_changed = 0;
for i = 1:N_chosen
    
    % Distance between chosen_X(i,:) and all solutions in pool
    [d, idx] = min(sqrt(sum((repmat(chosen_X(i,:), N_pool, 1) - pool.X) .^ 2, 2)));
    
    if d <= params.tol_ratio
        
        % Keep a copy to undo the replacement, if necessary
        bkp_x = pool.X(idx,:);
        bkp_y = pool.y(idx);
        bkp_age = pool.age(idx);
        
        % Replace the solution in the pool
        pool.X(idx,:) = chosen_X(i,:);
        pool.y(idx) = chosen_y(i);
        pool.age(idx) = max(pool.age) + 1;
        
        % Check the replacement
        if ~is_pool_ok(pool, params.tol_std)
            pool.X(idx,:) = bkp_x;
            pool.y(idx) = bkp_y;
            pool.age(idx) = bkp_age;
        else 
            d_pool_pop(idx) = min(sqrt(sum((repmat(pool.X(idx,:), N_pop, 1) - pop_X) .^ 2, 2)));
            has_pool_changed = 1;
        end
        
    else
        
        if N_pool < params.max_pool_size
            
            % Add the solution into the pool
            idx = N_pool + 1;
            pool.X(idx, :) = chosen_X(i,:);
            pool.y(idx) = chosen_y(i);
            pool.age(idx) = max(pool.age) + 1;
            
            % Check insertion
            if ~is_pool_ok(pool, params.tol_std)
                pool.X(idx,:) = [];
                pool.y(idx) = [];
                pool.age(idx) = [];
            else
                N_pool = N_pool + 1;
                d_pool_pop(idx) = min(sqrt(sum((repmat(pool.X(idx,:), N_pop, 1) - pop_X) .^ 2, 2)));
                has_pool_changed = 1;
            end
            
        else
            
            % Find the solution in the pool that is the farthest from pop_X
            [value, idx] = max(d_pool_pop);
            
            % Find the closest distance from chosen_X(i) to pop_X
            [d, idx_pop] = min(sqrt(sum((repmat(chosen_X(i,:), N_pop, 1) - pop_X) .^ 2, 2)));
            
            if d < value
                
                % Keep a copy to undo the replacement, if necessary
                bkp_x = pool.X(idx,:);
                bkp_y = pool.y(idx);
                bkp_age = pool.age(idx);

                % Replace the solution in the pool
                pool.X(idx,:) = chosen_X(i,:);
                pool.y(idx) = chosen_y(i);
                pool.age(idx) = max(pool.age) + 1;

                % Check the replacement
                if ~is_pool_ok(pool, params.tol_std)
                    pool.X(idx,:) = bkp_x;
                    pool.y(idx) = bkp_y;
                    pool.age(idx) = bkp_age;
                else 
                    d_pool_pop(idx) = d;
                    has_pool_changed = 1;
                end
            end
            
        end
        
    end
end

% Update age status in pool structure
if has_pool_changed
    pool.next_age = max(pool.age) + 1;
end

end
%
function [pool] = keeping_newest(chosen_X, chosen_y, pool, params)
% KEEPING_NEWEST: Upadate the pool of solutions.
% This function updates the pool by keeping the newest solutions evaluated
% in original function.

% Pool size N_pool
[N_pool, n] = size(pool.X);

% Remove repeated solutions from chosen_X
[chosen_X, idx, ~] = unique(chosen_X, 'rows', 'stable');
chosen_y = chosen_y(idx);
[N_chosen, ~] = size(chosen_X);

% Try to insert solutions in chosen_X into pool
has_pool_changed = 0;
for i = 1:N_chosen
    
    % Distance between chosen_X(i,:) and all solutions in pool
    [d, idx] = min(sqrt(sum((repmat(chosen_X(i,:), N_pool, 1) - pool.X) .^ 2, 2)));
    
    if d <= params.tol_ratio
        
        % Keep a copy to undo the replacement, if necessary
        bkp_x = pool.X(idx,:);
        bkp_y = pool.y(idx);
        bkp_age = pool.age(idx);
        
        % Replace the solution in the pool
        pool.X(idx,:) = chosen_X(i,:);
        pool.y(idx) = chosen_y(i);
        pool.age(idx) = max(pool.age) + 1;
        
        % Check the replacement
        if ~is_pool_ok(pool, params.tol_std)
            pool.X(idx,:) = bkp_x;
            pool.y(idx) = bkp_y;
            pool.age(idx) = bkp_age;
        else 
            has_pool_changed = 1;
        end
        
    else
        
        if N_pool < params.max_pool_size
            
            % Add the solution into the pool
            idx = N_pool + 1;
            pool.X(idx, :) = chosen_X(i,:);
            pool.y(idx) = chosen_y(i);
            pool.age(idx) = max(pool.age) + 1;
            
            % Check insertion
            if ~is_pool_ok(pool, params.tol_std)
                pool.X(idx,:) = [];
                pool.y(idx) = [];
                pool.age(idx) = [];
            else
                N_pool = N_pool + 1;
                has_pool_changed = 1;
            end
            
        else
            % Find the oldest solution in pool
            [value, idx] = min(pool.age);
            
            % Keep a copy to undo the replacement, if necessary
            bkp_x = pool.X(idx,:);
            bkp_y = pool.y(idx);
            bkp_age = pool.age(idx);

            % Replace the solution in the pool
            pool.X(idx,:) = chosen_X(i,:);
            pool.y(idx) = chosen_y(i);
            pool.age(idx) = max(pool.age) + 1;

            % Check the replacement
            if ~is_pool_ok(pool, params.tol_std)
                pool.X(idx,:) = bkp_x;
                pool.y(idx) = bkp_y;
                pool.age(idx) = bkp_age;
            else 
                has_pool_changed = 1;
            end
         
        end
        
    end
end

% Update age status in pool structure
if has_pool_changed
    pool.next_age = max(pool.age) + 1;
end

end
%
