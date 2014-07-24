function [zd,cd] = hp2disk(w,beta,z,c)
%HP2DISK Convert solution from the half-plane to one from the disk.
%   [Z,C] = HP2DISK(W,BETA,Z,C) quickly transforms the solution Z,C of
%   the Schwarz-Christoffel half-plane mapping parameter problem to the
%   solution ZD,CD of the disk problem.
%   
%   See also DISK2HP, HPPARAM, DPARAM.
 
%   Copyright 1998 by Toby Driscoll.
%   $Id: hp2disk.m 298 2009-09-15 14:36:37Z driscoll $

n = length(w);
zd = zeros(size(z));
if isinf(z(n))
  zd(n) = 1;
  zd(1:n-1) = (z(1:n-1)-i)./(z(1:n-1)+i);
else
  zd = (z-i)./(z+i);
  zd = zd/zd(n);
end
zd = sign(zd);

% Recalculate constant from scratch.
mid = (zd(1)+zd(2))/2;
qdat = scqdata(beta,16);
cd = (w(1) - w(2))/...
    (dquad(zd(2),mid,2,zd,beta,qdat) - dquad(zd(1),mid,1,zd,beta,qdat));

