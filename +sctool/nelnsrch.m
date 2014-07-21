function [retcode,xp,fp,Fp,maxtaken,nofun,btrack] = ...
            nelnsrch(xc,fc,fn,g,p,sx,sf,details,nofun,btrack,fparam)
%
%  [retcode,xp,fp,Fp,maxtaken,nofun,btrack] = ...
%         nelnsrch(xc,fc,fn,g,p,sx,sf,details,nofun,btrack,fparam)
%
%  This function is part of the Nonlinear Equations package and the
%  Unconstrained Minimization package, see NESOLVE.M or UMSOLVE.M.
%
%  It is a line search for use with Newton's Method of solving nonlinear
%  equations.  For function evaluations, it needs to know whether it is
%  doing Nonlinear Equations (NE) or Unconstrained Minimization (UM); it
%  distinguishes the two by the length of DETAILS, which is 16 for NE and
%  17 for UM.
%
%  Algorithm A6.3.1:  Part of the modular software system from the
%  appendix of the book "Numerical Methods for Unconstrained Optimization
%  and Nonlinear Equations" by Dennis & Schnabel, 1983.
%
%  Coded in MATLAB by Richard T. Behrens, March 1988.
%  Modified slightly for UM usage, January 1989.
%

%
% Initialization.
%
import sctool.*
n = length(xc);
xp = zeros(n,1);
fp = 0;
umflag = (length(details) == 17);     % This is how we tell NE from UM.

%
% Algorithm step 1.
%
maxtaken = 0;

%
% Algorithm step 2.
%
retcode = 2;

%
% Algorithm step 3.
%
alpha = 1E-4;

%
% Algorithm step 4.
%
newtlen = norm(sx .* p);

%
% Algorithm step 5.
%
if (newtlen > details(11))
   p = p * (details(11) / newtlen);
   newtlen = details(11);
end

%
% Algorithm step 6.
%
initslope = g'*p;

%
% Algorithm step 7.
%
rellength = max(abs(p)./max(abs(xc),(ones(n,1)./sx)));

%
% Algorithm step 8.
%
minlambda = details(9)/rellength;

%
% Algorithm step 9.
%
lambda = 1;

%
% Algorithm step 10.
%
bt = 0;
while (retcode >= 2)
   xp = xc + lambda*p;                      % step 10.1
   if umflag                                % step 10.2
      if details(15)
         fp = feval(fn,xp,fparam);
      else
         fp = feval(fn,xp);
      end
      nofun = nofun + 1;
   else
      if details(15)
         [fp,Fp,nofun] = nefn(xp,sf,fn,nofun,fparam);
      else
         [fp,Fp,nofun] = nefn(xp,sf,fn,nofun);
      end
   end
   if (fp <= fc + alpha*lambda*initslope)   % step 10.3a
      retcode = 0;
      maxtaken = ((lambda == 1) & (newtlen > 0.99*details(11)));
   elseif (lambda < minlambda)              % step 10.3b
      retcode = 1;
      xp = xc;
   else                                     % step 10.3c
      if (lambda == 1)
         if (details(1) > 0), disp('Quadratic Backtrack.'), end
         bt = bt + 1;
         lambdatemp = -initslope / (2*(fp-fc-initslope));
      else
         if (details(1) > 0), disp('Cubic Backtrack.'), end
         bt = bt + 1;
         a = (1/(lambda - lambdaprev)) * [1/lambda^2  (-1/lambdaprev^2); ...
               (-lambdaprev/(lambda^2))  lambda/(lambdaprev^2)] * ...
               [(fp-fc-lambda*initslope); (fpprev-fc-lambdaprev*initslope)];
         disc = a(2)^2 - 3*a(1)*initslope;
         if (a(1)==0)
            lambdatemp = -initslope/(2*a(2));
         else
            lambdatemp = (-a(2)+sqrt(disc))/(3*a(1));
         end
         if (lambdatemp > 0.5*lambda)
            lambdatemp = 0.5*lambda;
         end
      end
      lambdaprev = lambda;
      fpprev = fp;
      if (lambdatemp <= 0.1*lambda)
         lambda = 0.1*lambda;
      else
         lambda = lambdatemp;
      end
   end
end
if (bt < length(btrack))
   btrack(bt+1) = btrack(bt+1) + 1;
else
   btrack(bt+1) = 1;
end

