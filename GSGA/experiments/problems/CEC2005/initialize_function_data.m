function initialize_function_data
% Variables number
nvar = [10 20 30 50 100 150 200];
% Functions name
function_name = {'rastrigin_rot_func','hybrid_rot_func1','hybrid_rot_func2','hybrid_rot_func2_narrow',...
'sphere_func','rosenbrock_func','hybrid_composition','shifted_schwefel_problem_1_2'};

for f = 1:length(function_name)
    fname = function_name{f};
    save_function_data(fname,nvar)
end



end