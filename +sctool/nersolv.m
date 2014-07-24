function b = nersolv(M,M2,b)
%
%  b = nersolv(M,M2,b)
%
%  This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
%  It is a linear equation solve function for upper triangular systems.
%
%  Algorithm A3.2.2a:  Part of the modular software system from the
%  appendix of the book "Numerical Methods for Unconstrained Optimization
%  and Nonlinear Equations" by Dennis & Schnabel, 1983.
%
%  Coded in MATLAB by Richard T. Behrens, March 1988.
%

%
% Algorithm step 1.
%
n = length(M);
b(n) = b(n) / M2(n);

%
% Algorithm step 2.
%
for i = (n-1):-1:1
   b(i) = (b(i) - M(i,(i+1):n) * b((i+1):n)) / M2(i);
end

