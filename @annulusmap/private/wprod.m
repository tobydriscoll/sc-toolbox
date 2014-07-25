function w_prod = wprod(w,u,uw0,u_w1,dataz,param4)
% WPROD computes the product (DSC integrand)
% The definition of the integrand can be found in 
% USER'S GUIDE to DSCPACK, Chenglie Hu 1995

w = w(:);

mattemp = w(:,ones(length(uw0),1)) ./ uw0(ones(length(w),1),:);
param4temp = param4;
param4temp.POWER = param4temp.POWER(ones(length(w),1),ones(length(uw0),1),:);
param4temp.UARY = param4temp.UARY(ones(length(w),1),ones(length(uw0),1),:);
wth = log(wtheta(u,param4temp,mattemp));
w_prod = wth*(dataz.ALFA0(:)-1);

mattemp = w(:,ones(length(u_w1),1)) .* u_w1(ones(length(w),1),:);
param4temp = param4;
param4temp.POWER = param4temp.POWER(ones(length(w),1),ones(length(u_w1),1),:);
param4temp.UARY = param4temp.UARY(ones(length(w),1),ones(length(u_w1),1),:);
wth = log(wtheta(u,param4temp,mattemp));
w_prod = w_prod + wth*(dataz.ALFA1(:)-1);

w_prod = exp(w_prod);