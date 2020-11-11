function P = create_offsprings(pop_X, best_x, lb, ub, offsprings_per_solution, Nf, Nc, Ngood)
% CREATE_OFFSPRINGS: Create offspring solutions for each solution
% in the population.
%
% Input: 
%   pop_X: Current population of EA (rows are entries and coluns are the 
%       variables).
%   x_best: Best solution so far
%   lb: Lower bounds
%   ub: Upper bounds
%   offsprings_per_solution: Number of new solutions per solution
%   Nf:
%   Nc:
%   Ngood:
% Output:
%   P : Offsprings

switch Ngood
    case 0 
        if Nf < Nc
            [P] = create_offsprings_de_current_best1(pop_X, best_x, lb, ub, offsprings_per_solution);
        else
            [P] = create_offsprings_de_randbest(pop_X, best_x, lb, ub, offsprings_per_solution);
        end        
    otherwise
        [P] = create_offsprings_de_best(pop_X, best_x, lb, ub, offsprings_per_solution);        
end



end