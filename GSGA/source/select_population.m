function [pop_X,pop_y] = select_population(X,y,N)
% SELECT_POPULATION: Select solutions from the sample to
% compose the population of the SAEA.
%
% Input:
%   X: Sample (rows are entries and coluns are the variables)
%   y: Evaluation of each row in X
%   N: Population size.
%
% Output:
%   pop_X: Population of EA
%   pop_y: Evaluation of each row in pop_X
%

% Choose the N best solution 
[~,idx] = sort(y,'ascend');
pop_X = X(idx(1:N),:);
pop_y = y(idx(1:N));

end