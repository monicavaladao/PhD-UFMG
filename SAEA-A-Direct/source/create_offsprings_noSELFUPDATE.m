function [P, V, CR, Nf, DE_counter, any_DE_used] = create_offsprings_noSELFUPDATE(pop_X, pop_y, pop_V, pop_CR, best_x, lb, ub, ...
                                     DE_counter, bestever, offsprings_per_solution, Nf, Nc, Ngood, any_DE_used)
% CREATE_OFFSPRINGS: Create offspring solutions for each solution
% in the population.
%
% Input: 
%   pop_X: Current population of EA (rows are entries and coluns are the 
%       variables).
%   x_best: Best solution of the current population
%   lb: Lower bounds
%   ub: Upper bounds
%   offsprings_per_solution: Number of new solutions per solution
%   Nf: Feedback Mechanism parameters
%   Nc: Feedback Mechanism parameters
%   Ngood: Feedback Mechanism parameters
% Output:
%   P : Offsprings
%   V: Behavor correction of each ofsspring in P
%   Nf: Feedback Mechanism parameters
%   DE_counter: Counter of the number times in which the DE operator is used

% Number of solutions and variables 
[N,n] = size(pop_X);

% Define parameter
if n <= 100
    Nmod = 5;
else
    Nmod = 10;
end

V = pop_V;
CR = pop_CR;

switch Ngood
    case 0 
        if Nf < Nc
            [P, V, CR] = create_offsprings_slpso(pop_X, pop_y, pop_V, pop_CR, lb, ub, bestever);
            any_DE_used = 0;
        elseif (Nf >= Nc && Nf < 2*Nc)  
            if mod(Nf,Nmod) ~= 0
                pop_CR = 0.6*ones(N,1);
                [P] = create_offsprings_de_randbest(pop_X, best_x, lb, ub, pop_CR, offsprings_per_solution);
            else
                [P] = create_offsprings_de_current_best2(pop_X, best_x, lb, ub, pop_CR, offsprings_per_solution);
            end
            DE_counter = DE_counter + 1;
            any_DE_used = 1;
        else
            pop_CR = 0.8*ones(N,1);
            [P] = create_offsprings_de_rand(pop_X, lb, ub, pop_CR, offsprings_per_solution);
            DE_counter = DE_counter + 1;
            any_DE_used = 1;
            Nf = 1;            
       end        
    otherwise
        pop_CR = 0.2*ones(N,1);
        [P] = create_offsprings_de_current_best1(pop_X, best_x, lb, ub, pop_CR, offsprings_per_solution); 
        DE_counter = DE_counter + 1;
        any_DE_used = 1;
end

