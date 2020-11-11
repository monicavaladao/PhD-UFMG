function [exp_improv,y_hat] = expimp_SRGTSToolbox(x,y,info)
% EXPIMP_OODACE: Adapt "expimp" frunction from Alexander Forrester [1] to
% calculate the expected improvement using the metamodel built with SRGTSToolbox
%
% Input:
%   x:
%   y: 
%   info:
%
% Output:
%   exp_improv:
%   y_hat:
%
% References:
%   [1] FORESTER, Alexander. .....


% Best value in y so far
y_min = min(y);

% Prediction and MSE
[y_hat, predvar] = dace_predictor(x, info.SRTSToolbox.KRG_DACEModel);
mse = predvar;

% Expected Improvement
if mse == 0
    exp_improv = 0;
else
    ei_termone = (y_min - y_hat)*(0.5 + 0.5*erf((1/sqrt(2))*...
        ((y_min - y_hat)/sqrt(abs(mse)))));

    ei_termtwo = sqrt(abs(mse))*(1/sqrt(2*pi))*exp(-(1/2)*...
        ((y_min - y_hat)^2/mse));
    
    exp_improv = ei_termone + ei_termtwo;
end

end