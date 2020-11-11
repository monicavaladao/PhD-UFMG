function P = create_offsprings_de_randbest(pop_X, best_x, lb, ub, offsprings_per_solution)
% CREATE_OFFSPRINGS_DE_RANDBEST: Create offspring solutions for each solution
% in the population.
%
% Input: 
%   pop_X: Current population of EA (rows are entries and coluns are the 
%       variables).
%   x_best: Best solution so far
%   lb: Lower bounds
%   ub: Upper bounds
% Output:
%   P : Offsprings

% Population size (N) and number of variables (n)
[N, n] = size(pop_X);


% Generate new solutions
P = [];
for i = 1:N
    
    % Crossover rate
    C = 0.6;
    
    % Create for each solution
	U = zeros(offsprings_per_solution, n);
    F = struct();
    M = struct();
    
    for m = 1:offsprings_per_solution
        
        % Scale Factor
        %F = 0.8 + 0.2 * rand();
        %F = 0.4 + 0.5 * rand();
        F = 0.8;
        
        % Mutation strategy
        % DE/current-to-randbest/1:

        % Choose two solutions different from i
        j = randperm(N);
        while (sum(i == j(1:2)) > 0)
            j = randperm(N);
        end
        % Differential mutataion
        x_dif = pop_X(i,:) + rand* (best_x - pop_X(i,:)) + F*(pop_X(j(1),:) - pop_X(j(2),:));
        M(m).mutation_type = 'DE_best'; 
              
        % Truncation
        x_dif(x_dif < lb) = lb(x_dif < lb);
        x_dif(x_dif > ub) = ub(x_dif > ub);
        
        % Discrete recombination
        rand_idx = rand(n,1) <= C;
        rand_idx(randi(n)) = 1;
        U(m,:) = pop_X(i,:);
        U(m,rand_idx) = x_dif(1,rand_idx);		
    end
    
    P = [P;U];
    
end

end
