function [M,M1,M2,sing] = neqrdcmp(M)
%
%  [M,M1,M2,sing] = neqrdcmp(M)
%
%  This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
%  It is a QR decomposition function.  It differs from the one built
%  into MATLAB in that the result is encoded as rotation angles.  Also,
%  it is designed for square matrices only.
%
%  Algorithm A3.2.1:  Part of the modular software system from the
%  appendix of the book "Numerical Methods for Unconstrained Optimization
%  and Nonlinear Equations" by Dennis & Schnabel, 1983.
%
%  Coded in MATLAB by Richard T. Behrens, March 1988.
%

%
% Check size of input argument and allocate variables.
%
n = length(M);
M1 = zeros(n,1);
M2 = zeros(n,1);

%
% Algorithm step 1.
sing = 0;

%
% Algorithm step 2.
%
for k = 1:(n-1)
   eta = max(M(k:n,k));
   if (eta == 0)
      M1(k) = 0;
      M2(k) = 0;
      sing = 1;
   else
      M(k:n,k) = M(k:n,k) / eta;
      sigma = (sign(M(k,k))+(M(k,k)==0)) * norm(M(k:n,k));
      M(k,k) = M(k,k) + sigma;
      M1(k) = sigma * M(k,k);
      M2(k) = -eta * sigma;
      tau = (M(k:n,k)' * M(k:n,(k+1):n)) / M1(k);
      M(k:n,(k+1):n) = M(k:n,(k+1):n) - M(k:n,k) * tau;
   end
end

%
% Algorithm step 3.
%
M2(n) = M(n,n);

