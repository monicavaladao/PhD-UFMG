function [CR, CRm] = crossover_rate(CRm, N)
%  CROSSOVER_RATE: Create a crossover rate for each solution i = 1,...,N
% Input:
%   CRm: Parameter to update each CRi
%   N: Number of solutions
% Output:
%   CR: Vector crossover rate (Nx1)
%   CRm: Parameter to update each CRi

% CR
CR = normrnd(CRm,0.1,N,1);


end