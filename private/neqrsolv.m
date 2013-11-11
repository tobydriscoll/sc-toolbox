function b = neqrsolv(M,M1,M2,b)
%
%  b = neqrsolv(M,M1,M2,b)
%
%  This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
%  It is a linear equation solve function using the QR decomposition.
%
%  Algorithm A3.2.2:  Part of the modular software system from the
%  appendix of the book "Numerical Methods for Unconstrained Optimization
%  and Nonlinear Equations" by Dennis & Schnabel, 1983.
%
%  Coded in MATLAB by Richard T. Behrens, March 1988.
%

%
% Algorithm step 1.
%
n = length(M);
for j = 1:(n-1)
   tau = (M(j:n,j)' * b(j:n)) / M1(j);
   b(j:n) = b(j:n) - tau * M(j:n,j);
end

%
% Algorithm step 2.
%
b = nersolv(M,M2,b);

