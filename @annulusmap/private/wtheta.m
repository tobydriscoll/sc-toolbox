function w_theta = wtheta(u,param4,w)
% WTHETA evaluates theta-function at matrix w, 
% where u is the inner radius of the annulus. The 
% definition of theta function can be found in 
% USER'S GUIDE to DSCPACK, Chenglie Hu 1995
if (u < 0.63)   
    w_theta = -w;
    w_theta = w_theta(:,:,param4.ONEIU);
    w_theta = w_theta .^ param4.POWER; 
    w_theta = 1 + sum(param4.UARY .* (w_theta + 1./w_theta),3);
else
    wt = -i*log(-w);
    if (u >= 0.94)
        w_theta = exp(-0.25*(wt.^2)/(pi*param4.DLAM))/sqrt(param4.DLAM);
        return;
    end
    w_theta = 1 + 2*param4.VARY*cosh(wt/param4.DLAM);
    w_theta = exp(-0.25*(wt.^2)/(pi*param4.DLAM)).*(w_theta/sqrt(param4.DLAM));    
end