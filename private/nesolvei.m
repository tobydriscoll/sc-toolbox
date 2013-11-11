function [xf,termcode,path] = nesolvei(fvec,x0,details,fparam,jac,scale)
%FSOLVEI Solution to nonlinear equations with no initial Jacobian.
%	
%	Identical to FSOLVE, but the initial Jacobian is set to the identity
%	rather than being computed by finite differences. Obviously, you
%	need a special type of system to benefit from this change. The CRDT
%	system is one such.

%       Copyright (c) 1988 by the MathWorks, Inc.
%       Coded in MATLAB by Richard T. Behrens, April 1988.
%       Revised 11/27/88 JNL
%       Hookstep option added 8/21/90 RTB.

%	Modification by Toby Driscoll, 09/17/96.

% Initialization.
%
if (nargin < 3)
   details = zeros(16,1);
end
i = length(details);
if (i < 16)
   details((i+1):16,1) = zeros(16-i,1);  % For unspecified details.
end
if (nargin < 4)
   details(15) = 0;                      % No parameters to pass.
elseif isempty(fparam)
   details(15) = 0;
else
   details(15) = 1;
end
if (details(15) == 0)
   fparam = [];
end
if (nargin < 5)
   details(4) = 0;                       % No analytic jacobian.
elseif isempty(jac)
   details(4) = 0;
end
if (details(4) == 0)
   jac = '';
end
if (nargin < 6)
   scale = [];                           % No scaling given.
   if (details(16) == 2)
      details(16) = 1;
   end
end
if (nargout < 3)
   details(14) = 0;                      % No path output.
else
   details(14) = 1;
end
if details(14), path = x0.'; end
if (details(1) == 2), btrack = []; end
nofun = 0;              % Number of function evaluations.
trustvars = zeros(4,1); % variables for trust region methods.


%
% Algorithm step 2.
%
if (details(16) > 0)     % Might need F(x0) for scaling.
   if details(15)
      [fc,FVplus,nofun] = nefn(x0,ones(length(x0),1),fvec,nofun,fparam);
   else
      [fc,FVplus,nofun] = nefn(x0,ones(length(x0),1),fvec,nofun);
   end
else
   FVplus = zeros(length(x0),1);
end
[details,Sx,SF,termcode] = neinck(x0,FVplus,details,scale);

%
% Algorithm step 3.
%
if (termcode < 0)
   xf = x0;
   if details(14), path = [path;xf.']; end
   return
end

%
% Algorithm step 4.
%
itncount = 0;

%
% Algorithm step 5.
%
if details(15)
   [fc,FVplus,nofun] = nefn(x0,SF,fvec,nofun,fparam);
else
   [fc,FVplus,nofun] = nefn(x0,SF,fvec,nofun);
end

%
% Algorithm step 6.
%
%[termcode,consecmax] = nestop0(x0,FVplus,SF,details(8));
consecmax = 0;
if (max(SF .* abs(FVplus)) <= 1E-2 * details(8))
   termcode = 1;
else
   termcode = 0;
end

%
% Algorithm step 7.
%
if (termcode > 0)
   xf = x0;
   if (details(14)), path = [path;xf.']; end
else
   if details(4)
      if details(15)
         [Jc,addfun] = feval(jac,x0,fparam);
      else
         [Jc,addfun] = feval(jac,x0);
      end
      nofun = nofun + addfun;
  else
     %*** changed line *** 
     %[Jc,nofun] = nefdjac(fvec,FVplus,x0,Sx,details,nofun,fparam);
     Jc = eye(length(x0));
   end
   gc = Jc' * (FVplus .* (SF.^2));
   FVc = FVplus;
end

%
% Algorithm step 8.
%
xc = x0;

%
% Algorithm step 9.
%
restart = 1;
norest = 0;

%
% Algorithm step 10 (iteration).
%
%%if (details(1) > 0), clc, end
if details(1)
  %disp('ITN    F-COUNT    MAX(ABS(F))')
  wbh = waitbar(0,'Solution progress:');
  normFV0 = max(abs(FVc));
end
  
while (termcode == 0)
   if details(1)
      %disp(sprintf('%3d %8d %16g',itncount,nofun,max(abs(FVc))))
      fracdone = log(max(abs(FVc))/normFV0)/log(details(8)/normFV0);
      waitbar(fracdone)
   end
   itncount = itncount + 1;
   if (details(4) | details(3) | (1-details(5)))
      [M,Hc,sN] = nemodel(FVc,Jc,gc,SF,Sx,details(2));
   else
      error('Factored model not implemented.')
%     [] = nemodfac();
   end
   if (details(2) == 1)
      [retcode,xplus,fplus,FVplus,maxtaken,nofun,btrack] = ...
      nelnsrch(xc,fc,fvec,gc,sN,Sx,SF,details,nofun,btrack,fparam);
   elseif (details(2) == 2)
      [retcode,xplus,fplus,FVplus,maxtaken,details,trustvars,nofun] = ...
         nehook(xc,fc,fvec,gc,tril(M),Hc,sN,Sx,SF,details,itncount,trustvars,...
         nofun,fparam);
   else
      error('Dogleg not implemented.')
%     [] = nedogdrv();
   end
   if ((retcode ~= 1) | (restart) | (details(4)) | (details(3)))
      if (details(4))
         if details(15)
            [Jc,addfun] = feval(jac,xplus,fparam);
         else
            [Jc,addfun] = feval(jac,xplus);
         end
         nofun = nofun + addfun;
      elseif (details(3))
         [Jc,nofun] = nefdjac(fvec,FVplus,xplus,Sx,details,nofun,fparam);
      elseif (details(5))
         error('Factored secant method not implemented.')
%        [] = nebroyf();
      else
         Jc = nebroyuf(Jc,xc,xplus,FVc,FVplus,Sx,details(13));  % Broyden update
      end
      if (details(5))
         error('Gradient calculation for factored method not implemented.')
         % Calculate gc using QR factorization (see book).
      else
         gc = Jc' * (FVplus .* (SF.^2));
      end
      [consecmax,termcode] = nestop(xc,xplus,FVplus,fplus,gc,Sx,SF,retcode,...
                                details,itncount,maxtaken,consecmax);
   end
   if (((retcode == 1) | (termcode == 2)) & (1-restart) & ...
                                     (1-details(4)) & (1-details(3)))
      [Jc,nofun] = nefdjac(fvec,FVc,xc,Sx,details,nofun,fparam);
      gc = Jc' * (FVc .* (SF.^2));
      if ((details(2) == 2) | (details(2) == 3))
         details(7) = -1;
      end
      restart = 1;
      norest = norest + 1;
      if termcode==2, termcode = 0; end %***added by TAD
   else
      if (termcode > 0)
         xf = xplus;
         if (details(14)), path = [path;xf.']; end
      else
         restart = 0;
         if (details(14)), path = [path;xplus.']; end
      end
      xc = xplus;
      fc = fplus;
      FVc = FVplus;
%%      if (details(1) > 0)
%%         clc
%%         disp('The current iteration is: ')
%%         xc
%%      end
   end
end

if details(1)
  close(wbh)
end

if (details(1) == 2)
   %%disp('Press CR to see statistics . . .')
   %%fprintf([' ',7])
   %%pause

   %%clc
   %%format compact
   disp(' ')
   %%disp('Function: ')
   %%fvec
   %%disp('Starting point: ')
   %%x0.'
   %%disp('Termination condition: ')
   %%termcode
   disp(sprintf('Number of iterations: %d',itncount))
   disp(sprintf('Number of function evaluations: %d',nofun))
   disp(sprintf('Final norm(F(x)): %.6g',max(abs(FVc))))
   if ((1-details(3)) & (1-details(4)))
      disp(sprintf('Number of restarts for secant methods: %d',norest))
   end
%%   if (details(2) == 1)
%%      disp('Backtrack information: ')
%%      btrack
%%   end
%%   pause
end






