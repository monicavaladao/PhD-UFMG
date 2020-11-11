% Clear MATLAB workspace
clear all
close all
clc

% Remove folders from search path
rmpath('../source');
rmpath('../somtoolbox');
rmpath('./problems/CEC2010');
% Get problem data
%addpath('./problems/CEC2005');
%addpath('./problems/CEC2010');
%addpath('./problems/CEC2010/updated_datafiles');
%addpath('./problems/Other');
%addpath('./problems');
addpath(genpath('./problems'));
addpath('../source');
addpath('../somtoolbox');
%
addpath(genpath('../toolboxes/ooDACE'));
addpath(genpath('../toolboxes/SRGTSToolbox_DIRECT'));
addpath(genpath('../DIRECTalgorithm'));
%
%aux = load_analytic_problem('rastrigin_rot_func', 100);
%aux = load_analytic_problem('rastrigin_rot_func', 50);
%aux = load_analytic_problem('ackley', 100);
%aux = load_analytic_problem('ellipsoid', 100);
%aux = load_analytic_problem('ellipsoid', 50);
%aux = load_analytic_problem('rosen', 150);
%aux = load_analytic_problem('rosen', 200);
%aux = load_analytic_problem('ackley', 200);
%aux = load_analytic_problem('griewank', 100);
aux = load_analytic_problem('ackley', 50);
%aux = load_analytic_problem('rosenbrock_func', 100);
%aux = load_analytic_problem('rastrigin', 100);
%aux = load_analytic_problem('hybrid_rot_func1',50);
%aux = load_analytic_problem('rosenbrock_func', 100);
%aux = load_analytic_problem('rosenbrock_func', 20);
%aux = load_analytic_problem('sphere', 20);
%aux = load_analytic_problem('ackley', 10);
%aux = load_analytic_problem('rastrigin', 10);
%aux = load_analytic_problem('sphere_func', 100);
%aux = load_analytic_problem('griewank', 50);
%aux = load_analytic_problem('rosen', 50);
%aux = load_analytic_problem('rastrigin', 50);
%aux = load_analytic_problem('sg_shifted_rotated_rastrigin', 100);
%aux = load_analytic_problem('sg_shifted_rotated_elliptic', 50);
%aux = load_analytic_problem('sg_shifted_rotated_ackley', 20);
n = aux.n;
lb = aux.lb;
ub = aux.ub;
fobj = aux.fobj;

% Parameter of metamodel sample
rule = struct();
rule.name = 'newest';
rule.choose_sample = {'ChooseSample','newest', 'Verbose', false};

% Parameter of metamodel
metamodels = struct();

metamodels.name = 'ordinary-kriging';
metamodels.params = {'Metamodel', 'OrdinaryKriging_ooDACE', 'Verbose', false};

% metamodels.name = 'simple-kriging';
% metamodels.params = {'Metamodel', 'SimpleKriging_ooDACE', 'Verbose', false};

% 
% metamodels.name = 'universal-kriging1';
% metamodels.params = {'Metamodel', 'UniversalKriging1_ooDACE', 'Verbose', false};
% 
% metamodels.name = 'universal-kriging2';
% metamodels.params = {'Metamodel', 'UniversalKriging2_ooDACE', 'Verbose', false};
% 
% metamodels.name = 'blind-kriging';
% metamodels.params = {'Metamodel', 'BlindKriging_ooDACE', 'Verbose', false};

 %metamodels.name = 'rbf-gaussian';
 %metamodels.params = {'Metamodel', 'RBF_SRGTSToolbox', 'RBF', 'Gaussian', 'Verbose', false};
 %metamodels.name = 'rbf-multiquadric';
 %metamodels.params = {'Metamodel', 'RBF_SRGTSToolbox', 'RBF', 'Multiquadric', 'Verbose', false};


% Budget of objective function evaluations
max_eval = 2000;

% Size of initial sample
%npop = 50; % Size of each subpopulation
npop = 50 + floor(n/10); % Est� bom esse valor de npop
%npop = 80 + floor(n/10); % Est� bom esse valor de npop
%npop = 100 + floor(n/10);
ssize = 200;%100;%95; % Size of initial sample
X = lhsdesign(ssize, n);
X = repmat(lb, ssize, 1) + repmat(ub - lb, ssize, 1) .* X;
y = feval_all(fobj, X);

% Perform the TASEA algorithm 
[best_x, best_y, info] = surrogate_tasea_slpso(fobj, X, y, lb, ub, max_eval)
%[best_x, best_y, info] = surrogate_tasea_slpso_noDE(fobj, X, y, lb, ub, max_eval)
%[best_x, best_y, info] = surrogate_tasea_slpso_noSELFUPDATE(fobj, X, y, lb, ub, max_eval)


%[best_x, best_y, info] = surrogate_tasea_slpso(fobj, X, y, lb, ub, max_eval, metamodels.params{:}, rule.choose_sample{:})

% Remove folders from search path
rmpath('../source');
rmpath('../somtoolbox');
