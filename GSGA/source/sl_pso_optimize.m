function [xrbf_min,yrbf_min] = sl_pso_optimize(fobj,lb,ub)
%function [xrbf_min,yrbf_min] = sl_pso_optimize(fobj,lb,ub,numfobj)
% SL_PSO_OPTIMIZE: Adaptation of SL-PSO in order to optimize a RBF
% metamodel.
% Input:
%  fobj: Objective function (in this case it is the prediction function)
%  lb: Lower bound
%  up: Upper bound
%  numfobj: The maximal number of evaluations to be performed on fobj
% Output:
%   xrbf_min: Best solution
%   yrbf_min: Value of the best solution

% Initial particles
n = length(lb);
ssize = 50 + floor(n/10);
numfobj = 20*n;%10*n;%;
P = lhsdesign(ssize, n);
P = repmat(lb, ssize, 1) + repmat(ub - lb, ssize, 1) .*P;
FP = feval_all(fobj, P);

% Optimization process
[xrbf_min,yrbf_min,G,GXMIN,GOBJMIN,GNfobj,numfobj,RUNTime] = sl_pso(fobj,P,FP,lb,ub,ssize,numfobj);


end