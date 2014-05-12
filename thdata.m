function param4 = thdata(u)
%   THDATA generates data related only to inner radius u and used in
%   computing the theta-function.
param4 = struct('UARY','VARY','DLAM','IU','POWER','ONEIU');

if (u >= 0.63)
    param4.VARY = exp(pi^2/log(u));
    param4.DLAM = -log(u)/pi;
    return;
elseif (u < 0.06)
    param4.IU = 3;
elseif (u < 0.19)
    param4.IU = 4;
elseif (u < 0.33)
    param4.IU = 5;
elseif (u < 0.45)
    param4.IU = 6;
elseif (u < 0.55)
    param4.IU = 7;
else
    param4.IU = 8;
end

% Setting some parameter in advanced to be used in some functions.
param4.POWER = [];
param4.POWER(1,1,:) = (1:param4.IU);
param4.UARY = [];
param4.UARY(1,1,:) = (u*ones(1,param4.IU)).^((1:1:param4.IU).^2);
param4.ONEIU = ones(param4.IU,1);