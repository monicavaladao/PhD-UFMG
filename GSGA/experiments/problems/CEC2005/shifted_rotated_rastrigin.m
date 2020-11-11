% F10: Shifted Rotated Rastrigin's Function
% Properties:
%   Non-separable
%   Multimodal
%   Shifted
%   Scalable
%   Domain: [-5,5]^D
%   Local optima's number is huge
%   Global optmum: x = o(1:D), F10(x) = f_bias10 = -330
function [y] = shifted_rotated_rastrigin(x)
global initial_flag
initial_flag = 0;

func_num = 10;
fhd=str2func('rastrigin_rot_func'); %[-5,5]

load fbias_data;

y = feval(fhd,x)+f_bias(func_num);

% Turns the global optimum as zero
y = y - f_bias(func_num);

clear initial_flag
clear persistent 
end
%-----------------------------------------------
% 	10.Shifted Rotated Rastrign's Function 
function f=rastrigin_rot_func(x)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load rastrigin_func_data
    if length(o)>=D
        o=o(1:D);
        file_name = strcat('rastrigin_M','_D', int2str(D),'.mat');
        load(file_name);
    else
        c = 2;
        o=-5+10*rand(1,D);
        M=rot_matrix(D,c);
    end
 
    initial_flag=1;
end
x=x-repmat(o,ps,1);
x=x*M;
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end