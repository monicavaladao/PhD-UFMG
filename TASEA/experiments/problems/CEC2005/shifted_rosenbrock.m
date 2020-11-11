% F6: Shifted Rosenbrock's Function
% Properties:
%   Non-separable
%   Multimodal
%   Shifted
%   Scalable
%   Domain: [-100,100]^D
%   Having a very narrow valley from local optimum to global optimum
%   Global optmum: x = o(1:D), F6(x) = f_bias6 = 390
function [y] = shifted_rosenbrock(x)
global initial_flag
initial_flag = 0;

func_num = 6;
fhd=str2func('rosenbrock_func');

load fbias_data;

y = feval(fhd,x)+f_bias(func_num);

% Turns the global optimum as zero
y = y - f_bias(func_num);

clear initial_flag
clear persistent 
end
%--------------------------------------------------------
% 	6.Shifted Rosenbrock's Function
function f=rosenbrock_func(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load rosenbrock_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-90+180*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1)+1;
f=sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);
end