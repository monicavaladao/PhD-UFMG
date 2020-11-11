function [problem] = load_analytic_problem(name, n)
% LOAD_ANALYTIC_PROBLEM: load a problem by name and number of variables.
% 
% Input:
%	name: name of the problem
%	n   : number of variables
%
% Output:
%
% Problems available:
%   spheref: Sphere Function
%   sumpower: Sum of Differenr Power Function
% 	ellipsoid: Elipsoid Function
%	rosen: Rosen Function
%	ackley: Ackley Function
%	griewank: Griewank Function
%   rastrigin: Rastrigin Function
%	shifted_rotated_rastrigin: Shifted Rotated Rastrigin Function
%   rotated_hybrid_composition1: Rotated Version of Hybrid Composition Function
%   rotated_hybrid_composition_narrow: Rotated Hybrid Composition Function with narrow basin global
%   shifted_sphere: Shifted Sphere Function
%   shifted_rosenbrock: Shifted Rosenbrock Function


problem = struct();

switch name
    case 'sphere'
        problem.fobj = @spheref;
        problem.n = n;
        problem.name = 'Sphere';
        problem.lb = repmat(-5.12, 1, problem.n);
        problem.ub = repmat(5.12, 1, problem.n);
        
    case 'sumpower'
        problem.fobj = @sumpow;
        problem.n = n;
        problem.name = 'SumPower';
        problem.lb = repmat(-1, 1, problem.n);
        problem.ub = repmat(1, 1, problem.n);

   case 'ellipsoid'
		problem.fobj = @ellipsoid;
		problem.n = n;
		problem.name = 'Ellipsoid Function';
		problem.lb = repmat(-5.12, 1, problem.n);
		problem.ub = repmat(5.12, 1, problem.n); 
	
    case 'rosen'
		problem.fobj = @rosen;
		problem.n = n;
		problem.name = 'Rosen Function';
		problem.lb = repmat(-2.048, 1, problem.n);
		problem.ub = repmat(2.048, 1, problem.n);

 case 'ackley'
		problem.fobj =@ackley;
		problem.n = n;
		problem.name = 'Ackley Function';
		problem.lb = repmat(-32.768, 1, problem.n);
		problem.ub = repmat(32.768, 1, problem.n);   

	case 'griewank'
		problem.fobj = @griewank;
		problem.n = n;
		problem.name = 'Griewank Function';
		problem.lb = repmat(-600, 1, problem.n);
		problem.ub = repmat(600, 1, problem.n);	
		       
     case 'rastrigin'
		problem.fobj = @rastrigin;
		problem.n = n;
		problem.name = 'Rastrigin';
		problem.lb = repmat(-5.12, 1, problem.n);
		problem.ub = repmat(5.12, 1, problem.n);  
        
     %----------------CEC2005---------------------      
     case 'rastrigin_rot_func' % Shifted Rotated Rastrigin Function
		problem.fobj = @shifted_rotated_rastrigin;
		problem.n = n;
		problem.name = 'ShiRotRast';
		problem.lb = repmat(-5, 1, problem.n);
		problem.ub = repmat(5, 1, problem.n);

     case 'hybrid_rot_func1' % Rotated Version of Hybrid Composition Function
		problem.fobj = @rotated_hybrid_composition1;
		problem.n = n;
		problem.name = 'RotHybComp1';
		problem.lb = repmat(-5, 1, problem.n);
		problem.ub = repmat(5, 1, problem.n);
        
     case 'hybrid_rot_func2_narrow' % Rotated Hybrid Composition Function with narrow basin global
		problem.fobj = @rotated_hybrid_composition_narrow;
		problem.n = n;
		problem.name = 'RotHybComp2Narrow';
		problem.lb = repmat(-5, 1, problem.n);
		problem.ub = repmat(5, 1, problem.n);
        
      case 'sphere_func' % Shifted Sphere
		problem.fobj = @shifted_sphere;
		problem.n = n;
		problem.name = 'ShifSphere';
		problem.lb = repmat(-100, 1, problem.n);
		problem.ub = repmat(100, 1, problem.n);
        
     case 'rosenbrock_func' % Shifted Rosenbrock Function
		problem.fobj = @shifted_rosenbrock;
		problem.n = n;
		problem.name = 'ShifRosen';
		problem.lb = repmat(-100, 1, problem.n);
		problem.ub = repmat(100, 1, problem.n);

end