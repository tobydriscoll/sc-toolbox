function [consecmax,termcode] = nestop(xc,xp,F,Fnorm,g,sx,sf,retcode,...
                                details,itncount,maxtaken,consecmax)
%
%  [consecmax,termcode] = nestop(xc,xp,F,Fnorm,g,sx,sf,retcode,...
%                                details,itncount,maxtaken,consecmax)
%
%  This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
%  It decides whether or not to stop iterating when solving a set of
%  nonlinear equations.
%
%  Algorithm A7.2.3:  Part of the modular software system from the
%  appendix of the book "Numerical Methods for Unconstrained Optimization
%  and Nonlinear Equations" by Dennis & Schnabel, 1983.
%
%  Coded in MATLAB by Richard T. Behrens, March 1988.
%

%
% Algorithm step 1.
%
n = length(xc);
termcode = 0;

%
% Algorithm step 2.
%
if (retcode == 1)
   termcode = 3;
elseif (max(sf .* abs(F)) <= details(8))
   termcode = 1;
elseif (max( abs(xp - xc) ./ max(abs(xp),ones(n,1)./sx)) <= details(9))
   termcode = 2;
elseif (itncount >= details(6))
   termcode = 4;
elseif (maxtaken)
   consecmax = consecmax + 1;
   if (consecmax == 5)
      termcode = 5;
   end
else
   consecmax = 0;
   if (details(4) | details(3))
      if (max(abs(g).*max(abs(xp),ones(n,1)./sx)/max(Fnorm,(n/2))) ...
                                         <= details(10))
         termcode = 6;
      end
   end
end

