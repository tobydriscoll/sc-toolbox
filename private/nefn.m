function [fplus,FVplus,nofun] = nefn(xplus,SF,fvec,nofun,fparam)
%
%  [fplus,FVplus] = nefn(xplus,SF,fvec)
%
%  This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
%  It evaluates the vector function and calculates the sum of squares
%  for nonlinear equations.
%
%  Part of the modular software system from the appendix of the book
%  "Numerical Methods for Unconstrained Optimization and Nonlinear
%  Equations" by Dennis & Schnabel, 1983.
%
%  Coded in MATLAB by Richard T. Behrens, April 1988.
%

if (nargin < 5)
   FVplus = feval(fvec,xplus);
else
   FVplus = feval(fvec,xplus,fparam);
end
fplus = .5 * sum((SF .* FVplus).^2);
nofun = nofun + 1;

