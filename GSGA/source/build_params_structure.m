function [problem, params] = build_params_structure(fobj,X,y,lb,ub,max_eval,varargin)
% BUILD_PARAMS_STRUCTURE: Build an auxiliar structure used to keep 
% parameters used by the SAEA algorithm.
%
% Input:
%   fobj: handle to the objective function
%   X: Sample (rows are entries and coluns are the variables)
%   y: Evaluation of each row in X
%   lb: Lower bounds
%   ub: Upper bounds
%
% Optional input (key/value pairs):
%   - EvolutionControl: Strategy used to control solutions.
%       Values: 'metamodel', 'random'.
%   - Metamodel: Type of metamodel.
%       Values: 'OrdinaryKriging', 'UniversalKriging1',
%       'UniversalKriging2', 'BlindKriging', 'RBF'.
%   - Optimizer: Algorithm used to optimize the metamodel parameters. It is
%       used with Kriging metamodels only.
%       Values: 'sqp', 'fmincon', 'ga'.
%   - RBF: Type of RBF function. It is used with RBF metamodel only.
%       Values: 'Gaussian', 'GaussianCrossValidation', 'Multiquadric'.
%
% Output:
%   params: A structure with fields that keep parameters used by the SAEA
%       (including toolbox ooDACE).


% Sample size and number of variables
[init_sample_size, n] = size(X);

% Problem data
problem = struct();
problem.fobj = fobj;
problem.n = n;
problem.lb = lb;
problem.ub = ub;

% Parse optional input parameters
default_optimizer = 'fmincon';

parser = inputParser;
parser.addParameter('Verbose', true, @islogical);
%parser.addParameter('Metamodel_Kriging', 'OrdinaryKriging_ooDACE', @ischar);
%parser.addParameter('Metamodel_Kriging', 'UniversalKriging1_ooDACE', @ischar);
%parser.addParameter('Metamodel_Kriging', 'OrdinaryKriging_SRGTSToolbox', @ischar);
parser.addParameter('Metamodel_Kriging', 'UniversalKriging1_SRGTSToolbox', @ischar);
parser.addParameter('Metamodel_RBF', 'RBF_SRGTSToolbox', @ischar);
parser.addParameter('Optimizer', default_optimizer, @ischar);
%parser.addParameter('RBF', 'Gaussian', @ischar);
parser.addParameter('RBF', 'Multiquadric', @ischar);
parser.parse(varargin{:});

% Build parameters structure
params = struct();

% Budget of objective function evaluations
params.max_eval = max_eval;

% Parameters from parser
params.verbose = parser.Results.Verbose;
params.metamodel_kriging = parser.Results.Metamodel_Kriging;
params.metamodel_rbf = parser.Results.Metamodel_RBF;
params.optimizer = parser.Results.Optimizer;
params.rbf = parser.Results.RBF;

% Threshold values
params.tol_epsilon = 1e-1;
params.tol_ratio = 1e-6;
params.tol_std = 1e-6;

% Dimension of the reduced low-dimensional space
params.LDS.l = 4;

% Number of evaluations already spent
params.init_sample_size = init_sample_size;

% Set EA population size
params.pop_size = 50;

% Way to choose de metamodel sample
params.choose_sample.rbf = 'nearest';
params.choose_sample.kriging = 'newest';

% To updtade the pool of solutions
params.choose_sample.rule = 'newest';

% Set Metamodel sample size
params.sample_size.rbf = min(5*n + 10,init_sample_size);%200;%400;%100;
params.sample_size.kriging = init_sample_size;%400;%100;
params.sample_max_sample_size = 400;%500;%400;%200;%120; % The same as max_pool_size


if init_sample_size < params.sample_size.rbf
    error(strcat('Initial sample should have size(X,1)>=', ...
        int2str(params.sample_size), ' to be able to build the metamodel.'));
end

if init_sample_size < params.pop_size
    error(strcat('Initial sample should have size(X,1)>=', ...
        int2str(params.pop_size), ' to be able to select the population.'));
end

% Set the maximum number of solution in the pool of solutions evaluated
% with the original objective function
params.max_pool_size = max(init_sample_size,400);%400;%2*params.sample_size;

% ooDACE parameters
params.oodace = struct();

% MATLAB fmincon optimizer (used on Kriging metamodel)
nvar_opt = 1;
params.oodace.GradObj = 'on';
params.oodace.DerivativeCheck = 'off';
params.oodace.Diagnostics = 'off';
params.oodace.Algorithm = 'active-set';
params.oodace.MaxFunEvals = 100*nvar_opt;
params.oodace.MaxIter = 400;%100;
params.oodace.TolFun = 1e-4;
params.oodace.GradObj = 'off';
params.oodace.hpOptimizer = MatlabOptimizer(nvar_opt, 1, params.oodace);

% SRGTSToolbox parameters
params.srgts = struct();
% RBF Metamodel with Gaussian base function with sigma=1
params.srgts.P = X;
params.srgts.T = y;        
params.srgts.FIT_Fn = @rbf_build;
params.srgts.RBF_type = 'G';
params.srgts.RBF_c = 1;
params.srgts.RBF_usePolyPart = 0;
params.srgts.SRGT = 'RBF';


% GSGA paramenters
params.gsga.Ns = 10; % Number of selected solutions by prescreening based on Simplified Kriging
params.gsga.kmax = 3; 
params.gsga.N1 = round(params.pop_size*0.9); % Number of individuals that are created for crossover
if mod(params.gsga.N1,2) ~= 0
    params.gsga.N1 = params.gsga.N1 - 1;
end
params.gsga.N2 = round(params.pop_size*0.1);
params.gsga.N3 = params.pop_size - params.gsga.Ns; % Number of survided individuals
params.gsga.pcross = 0.9; % Crossover probability
params.gsga.pmut = 1; % Mutation probability


end