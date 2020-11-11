% F13: Shifted Expanded Griewan's plus Rosenbrock's Function (F8F2)
% Properties:
%   Non-separable
%   Multimodal
%   Shifted
%   Scalable
%   Domain: [-3,1]^D
%   Having a very narrow valley from local optimum to global optimum
%   Global optmum: x = o, F13(o) = f_bias13 = -130
function [y] = exp_ext_griewank_rosenbrock(x)
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
%-------------------------------------------------------
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
%--------------------------------------------------------
% 	13. Expanded Extended Griewank's plus Rosenbrock's Function (F8F2)
function fit=EF8F2_func(x)
%-3,1
global initial_flag
persistent  o 
[ps,D]=size(x);
if initial_flag==0
    load EF8F2_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-1+1*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1)+1;
fit=0;
for i=1:(D-1)
    fit=fit+F8F2(x(:,[i,i+1]));
end
    fit=fit+F8F2(x(:,[D,1]));
end
%---------------------------------------------------------------------    
function f=F8F2(x)
f2=100.*(x(:,1).^2-x(:,2)).^2+(1-x(:,1)).^2;
f=1+f2.^2./4000-cos(f2);
end