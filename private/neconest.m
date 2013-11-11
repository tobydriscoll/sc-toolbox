function est = neconest(M,M2)
%
%  This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
%  est = neconest(M,M2)
%  This is an estimate of the l-1 condition number of an upper triangular
%  matrix.
%
%  Algorithm A3.3.1:  Part of the modular software system from the
%  appendix of the book "Numerical Methods for Unconstrained Optimization
%  and Nonlinear Equations" by Dennis & Schnabel, 1983.
%
%  Coded in MATLAB by Richard T. Behrens, March 1988.
%

%
% Allocate variables.
%
n = length(M);
p = zeros(n,1);
pm = zeros(n,1);
x = zeros(n,1);

%
% Algorithm steps 1 & 2.
%
est = norm( triu(M)-diag(diag(M))+diag(M2) ,1);

%
% Algorithm step 3.
%
x(1) = 1/M2(1);

%
% Algorithm step 4.
%
p(2:n) = M(1,2:n) * x(1);

%
% Algorithm step 5.
%
for j = 2:n
   xp = (+1-p(j)) / M2(j);
   xm = (-1-p(j)) / M2(j);
   temp  = abs(xp);
   tempm = abs(xm);
   for i = (j+1):n
      pm(i) = p(i) + M(j,i)*xm;
      tempm = tempm + abs(p(i))/abs(M2(i));
      p(i) = p(i) + M(j,i) * xp;
      temp = temp + abs(p(i))/abs(M2(i));
   end
   if (temp > tempm)
      x(j) = xp;
   else
      x(j) = xm;
      p((j+1):n) = pm((j+1):n);
   end
end

%
% Algorithm steps 6 & 7.
%
est = est / norm(x,1);

%
% Algorithm step 8.
%
x = nersolv(M,M2,x);

%
% Algorithm steps 9 & 10.
%
est = est * norm(x,1);

