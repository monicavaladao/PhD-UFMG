% F2: Shifted Schwefel's Problem 1.2
% Properties:
%   Non-separable
%   Unimodal
%   Shifted
%   Scalable
%   Domain: [-100,100]^D
%   Global optmum: x = o, F2(o) = f_bias2 = -450

function [y] = shifted_schwefel_problem_1_2(x)
global initial_flag
initial_flag = 0;

func_num = 2;
fhd=str2func('schwefel_102'); %[-100,100]

load fbias_data;

y = feval(fhd,x)+f_bias(func_num);

% Turns the global optmum as zero
y = y - f_bias(func_num);

clear initial_flag
clear persistent 
end
%--------------------------------------------------------------------------
% 	2.Shifted Schwefel's Problem 1.2
function f=schwefel_102(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load schwefel_102_data
   if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1);
f=0;
for i=1:D
    f=f+sum(x(:,1:i),2).^2;
end
end
