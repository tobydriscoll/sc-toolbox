function [xp,fp,Fp,maxtaken,retcode,xpprev,fpprev,Fpprev,details,nofun] = ...
   netrust(retcode,xpprev,fpprev,Fpprev,xc,fc,fn,g,L,s,sx,sf,...
   newttaken,details,steptype,H,umflag,nofun,fparam)
%
% function [xp,fp,Fp,maxtaken,retcode,xpprev,fpprev,Fpprev,details,nofun] = ...
%          netrust(retcode,xpprev,fpprev,Fpprev,xc,fc,fn,g,L,s,sx,sf,...
%          newttaken,details,steptype,H,umflag,nofun,fparam)
%
%  This function is part of the Nonlinear Equations package and the
%  Unconstrained Minimization package, see NESOLVE.M or UMSOLVE.M.
%
%  It decides whether or not the proposed step is acceptable, and adjusts
%  the trust radius accordingly.
%

%
%  Algorithm A6.4.5:  Part of the modular software system from the
%  appendix of the book "Numerical Methods for Unconstrained Optimization
%  and Nonlinear Equations" by Dennis & Schnabel, 1983.
%
%  Coded in MATLAB by Richard T. Behrens, August 1990.
%

%
% Initialization.
%
n = length(xc);
xp = zeros(n,1);
fp = 0;

%
% Algorithm steps 1-4.
%
maxtaken = 0;
alpha = 1e-4;
steplen = norm(sx.*s);
xp = xc + s;

%
% Algorithm step 5.
%
if umflag
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

%
% Algorithm steps 6-8.
%
deltaf = fp - fc;
initslope = g'*s;
if (retcode~=3), fpprev = 0; end

%
% Algorithm step 9.
%
if ((retcode==3)&((fp>=fpprev)|(deltaf>alpha*initslope)))   % step 9a.
   retcode = 0;
   xp = xpprev;
   fp = fpprev;
   Fp = Fpprev;
   details(7) = details(7)/2;
%   if (details(1)>0), disp('Decreasing trust radius.'), end
else
   if (deltaf >= alpha*initslope)         % step 9b.
      rellength = max(abs(s)./max([abs(xp)'; 1.0./sx'])');
      if (rellength<details(9))       % details(9) is steptol
         retcode = 1;
         xp = xc;
      else
         retcode = 2;
         deltatemp = -initslope*steplen/(2*(deltaf-initslope));
         if (deltatemp<.1*details(7))
            details(7) = .1*details(7);
         elseif (deltatemp>.5*details(7))
            details(7) = .5*details(7);
         else
            details(7) = deltatemp;
         end
%         if (details(1)>0), disp('Decreasing trust radius.'), end
      end
   else                                   % step 9c.
      deltafpred = initslope;
      if (steptype==1)
         deltafpred = deltafpred + .5*(s'*H*s);
      else
         ttemp = L'*s;
         deltafpred = deltafpred + .5*(ttemp'*ttemp);
      end
      if ((retcode~=2)&((abs(deltafpred-deltaf)<=.1*abs(deltaf)) | ...
           (deltaf<=initslope)) & (~newttaken) & (details(7)<=.99*details(11)))
         retcode = 3;
         xpprev = xp;
         fpprev = fp;
         Fpprev = Fp;
         details(7) = min(2*details(7),details(11));  % details(11) is maxstep
%         if (details(1)>0), disp('Increasing trust radius.'), end
      else
         retcode = 0;
         if (steplen > .99*details(11)), maxtaken=1; end
         if (deltaf>=.1*deltafpred)
            details(7) = details(7)/2;
%            if (details(1)>0), disp('Decreasing trust radius.'), end
         elseif (deltaf<=.75*deltafpred)
            details(7) = min(2*details(7),details(11));
%            if (details(1)>0), disp('Increasing trust radius.'), end
         end
      end
   end
end

