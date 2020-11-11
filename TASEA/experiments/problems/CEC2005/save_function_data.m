function save_function_data(fname,nvar)
% Update some CEC05's functions according to the number of variables D
% Imput: 
%   fname: Name of the objetive function
%   nvar: Vector whith the number of variables
% Example:
%   save_function_data('rastrigin_rot_func',1000)


nD = length(nvar);
for i = 1:nD
    D = nvar(i);
    
    switch fname
        case 'rastrigin_rot_func'
            o = -5+10*rand(1,D);
            save('rastrigin_func_data','o');
            c = 2;
            M = rot_matrix(D,c);
            filename = strcat('rastrigin_M','_D', int2str(D),'.mat');
            save(filename,'M');
            %
        case 'hybrid_rot_func1'
            fun_num=10;
            o = -5+10*rand(fun_num,D);
            save('hybrid_func1_data.mat','o');
            c=[2,2,2,2,2,2,2,2,2,2];
            for i=1:fun_num
                eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
            end
            filename = strcat('hybrid_func1_M','_D', int2str(D),'.mat');
            save(filename,'M');
        case 'hybrid_rot_func2_narrow'
            fun_num=10;
            o=-5+10*rand(fun_num,D);
            save('hybrid_func2_data.mat','o');
            c=[2 3 2 3 2 3 20 30 200 300];
            for i=1:fun_num
                eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
            end
            filename = strcat('hybrid_func2_M','_D', int2str(D),'.mat');
            save(filename,'M');
        case 'hybrid_rot_func2'
            fun_num=10;
            o=-5+10*rand(fun_num,D);
            save('hybrid_func2_data.mat','o');
            c=[2 3 2 3 2 3 20 30 200 300];
            for i=1:fun_num
                eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
            end
            filename = strcat('hybrid_func2_M','_D', int2str(D),'.mat');
            save(filename,'M');
        case 'sphere_func'
            o=-100+200*rand(1,D);
            save('sphere_func_data.mat','o');
        case 'rosenbrock_func'
            o=-90+180*rand(1,D);
            save('rosenbrock_func_data.mat','o');
        case 'hybrid_composition'
            fun_num = 10;
            o = -5+10*rand(fun_num,D);
            save('hybrid_func1_data.mat','o');
       case 'shifted_schwefel_problem_1_2'
            o = -100+200*rand(1,D);
            save('schwefel_102_data.mat','o');   
            
            
    end
   
end


end