function [J,nofun] = nefdjac(fvec,fc,xc,sx,details,nofun,fparam)
%
% This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
% [J,nofun] = nefdjac(fvec,fc,xc,sx,details,nofun,fparam)
% This is a "Finite Differance Jacobian Approximation". It
% calculates a finite differance appproximation to J(xc)
% (the Jacobin of F(x) at x = xc).
%
% Algorithm A5.4.1: Part of the modular software system from
% the appendix of the book "Numerical Methods for Unconstrained
% Optimization and Nonlinear Equations" by Dennis & Schnabel 1983.
%
% Coded in Matlab by Sherkat Masoum M., March 1988.
% Edited by Richard T. Behrens, June 1988.
%

%
% Algorithm step 1.
%
n=length(fc);
sqrteta = sqrt(details(13));

%
% Algorithm step 2.
%
for j =1:n
  stepsizej = sqrteta * max(abs(xc(j)),1/sx(j)) * (sign(xc(j))+(xc(j)==0));
%   To incorporate a different stepsize rule, change the previous line.
  tempj = xc(j);
  xc(j) = xc(j) + stepsizej;
  stepsizej=xc(j)-tempj;
%   The previous line reduces finite precision error slightly,
%   see section 5.4 of the book.
  if details(15)
    fj =feval(fvec,xc,fparam);       % Evaluate function w/parameters.
  else
    fj =feval(fvec,xc);              % Evaluate function w/o parameters.
  end
  nofun = nofun + 1;
  J(1:n,j) = (fj(1:n) - fc(1:n))/stepsizej;
  xc(j) = tempj;
end

