% F1: Shifted Rotated Rastrigin's Function
% Properties:
%   Separable
%   Unimodal
%   Shifted
%   Scalable
%   Domain: [-100,100]^D
%   Global optimum: x = o(1:D), F1(x) = f_bias1 = -450
function [y] = shifted_sphere(x)
global initial_flag
initial_flag = 0;

func_num = 1;
fhd=str2func('sphere_func');

load fbias_data;

y = feval(fhd,x)+f_bias(func_num);

% Turns the global optimum as zero
y = y - f_bias(func_num);

clear initial_flag
clear persistent 
end
%----------------------------------------------------
% 	1.Shifted Sphere Function 
function fit=sphere_func(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load sphere_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1);
fit=sum(x.^2,2);

end