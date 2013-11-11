function [retcode,xp,fp,Fp,maxtaken,details,trustvars,nofun] = ...
          nehook(xc,fc,fn,g,L,H,sN,sx,sf,details,itn,trustvars,nofun,fparam)
%
% [retcode,xp,fp,Fp,maxtaken,details,trustvars,nofun] = ...
%    nehook(xc,fc,fn,g,L,H,sN,sx,sf,details,trustvars,nofun,fparam)
%
%  This function is part of the Nonlinear Equations package and the
%  Unconstrained Minimization package, see NESOLVE.M or UMSOLVE.M.
%
%  It is a driver for locally constrained optimal ("hook") steps for use
%  with Newton's Method of solving nonlinear equations.  For function
%  evaluations, it needs to know whether it is doing Nonlinear Equations
%  (NE) or Unconstrained Minimization (UM); it distinguishes the two by the
%  length of DETAILS, which is 16 for NE and 17 for UM.
%
%  TRUSTVARS is a vector of variables that, though not used externally,
%  need to be preserved between calls to NEHOOK.  The elements are
%  defined as:
%           1 = mu
%           2 = deltaprev
%           3 = phi
%           4 = phiprime

%
%  Algorithms A6.4.1 and A6.4.2:   Incorporates both the "hookdriver" and
%  "hookstep" algorithms.  Part of the modular software system from the
%  appendix of the book "Numerical Methods for Unconstrained Optimization
%  and Nonlinear Equations" by Dennis & Schnabel, 1983.
%
%  Coded in MATLAB by Richard T. Behrens, August 1990.
%

%
% Initialization.
%
n = length(xc);
umflag = (length(details) == 17);     % This is how we tell NE from UM.
xpprev = zeros(n,1);                  % allocation
fpprev = 0;                           % allocation
Fpprev = zeros(n,1);                  % allocation
%
% Algorithm steps 1-3.
%
retcode = 4;
firsthook = 1;
newtlen = norm(sx.*sN);

%
% Algorithm step 4.
%
if ((itn==1)|(details(7)==-1))        % details(7) is delta.
   trustvars(1) = 0;                  % trustvars(1) is mu.
   if (details(7)==-1)
      alpha = g./sx; alpha = alpha'*alpha;
      beta = L'*(g./(sx.*sx)); beta = beta'*beta;
      details(7) = (alpha^1.5)/beta;
      if (details(7) > details(11))   % details(11) is maxstep.
         details(7) = details(11);
      end
   end
end

%
% Algorithm step 5 (incorporating algorithm A6.4.2).
%
while (retcode >= 2)             % Calculate and check a new step.
   hi = 1.5; lo = 0.75;          % Start of A6.4.2.
   if (newtlen <= hi*details(7))
      newttaken = 1;
      s = sN;
      trustvars(1) = 0;          % trustvars(1) is mu.
      details(7) = min(details(7),newtlen);
   else
      newttaken = 0;
      if (trustvars(1) > 0)
         trustvars(1) = trustvars(1) - ...
            ((trustvars(3) + trustvars(2))/details(7))* ...
            (((trustvars(2)-details(7))+trustvars(3))/trustvars(4));
      end
      trustvars(3) = newtlen - details(7);
      if firsthook
         firsthook = 0;
         tempvec = L\((sx.*sx).*sN);
         phiprimeinit = -(tempvec'*tempvec)/newtlen;
      end
      mulow = -trustvars(3)/phiprimeinit;
      muup = norm(g./sx)/details(7);
      done = 0;
      while (~done)
         if ((trustvars(1) < mulow)|(trustvars(1)>muup))
            if (mulow<0), disp('warning, mulow<0'), keyboard, end
            trustvars(1) = max(sqrt(mulow*muup),muup*1e-3);
         end
         [L642,maxadd] = nechdcmp(H+trustvars(1)*diag(sx.*sx),0);
         s = -L642'\(L642\g);    % L642 is a copy of L local to A6.4.2.
         steplen = norm(sx.*s);
         trustvars(3) = steplen - details(7);
         tempvec = L642\((sx.*sx).*s);
         trustvars(4) = -(tempvec'*tempvec)/steplen;
         if (((steplen>=lo*details(7))&(steplen<=hi*details(7))) ...
                 | (muup-mulow<=0))
            done = 1;
         else
            mulow = max(mulow,trustvars(1)-(trustvars(3)/trustvars(4)));
            if (trustvars(3)<0), muup = trustvars(1); end
            trustvars(1) = trustvars(1) - ((steplen/details(7))* ...
                           (trustvars(3)/trustvars(4)));
         end
      end
   end                             % End of A6.4.2.
   trustvars(2) = details(7);      % trustvars(2) is deltaprev.
   [xp,fp,Fp,maxtaken,retcode,xpprev,fpprev,Fpprev,details,nofun] = ...
      netrust(retcode,xpprev,fpprev,Fpprev,xc,fc,fn,g,L,s,sx,sf,...
      newttaken,details,1,H,umflag,nofun,fparam);
end

