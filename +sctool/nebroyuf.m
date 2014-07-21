function A = nebroyuf(A,xc,xp,fc,fp,sx,eta)
%
% A = nebroyuf(A,xc,xf,fc,fp,sx,eta)
%
% This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
% This updates A, a secant approximation to the jacobian, using
% BROYDEN'S UNFACTORED SECANT UPDATE.
%
% Algorithm A8.3.1:  Part of the modular software system from
% the appendix of the book "Numerical Methods for Unconstrained
% Optimization and Nonlinear Equations" by Dennis & Schnabel 1983.
%
%
% Coded in Matlab by Sherkat Masoum M., April 1988.
% Edited by Richard T. Behrens, June 1988.
%

%
% Algorithm step 1.
%
n = length(A);
s=xp-xc;

%
% Algorithm step 2.
%
denom=norm(sx.*s)^2;

%
% Algorithm step 3.
%
tempi = (fp - fc - A*s);
ii = find(abs(tempi) < eta*(abs(fp)+abs(fc)));
tempi(ii) = zeros(length(ii),1);
A = A + (tempi/denom)*(s.*(sx.*sx))';

