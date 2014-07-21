function [dout,Sx,SF,termcode] = neinck(x0,F0,din,scale)
%
% [dout,Sx,SF,termcode] = neinck(x0,din,scale)
%
%  This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
%  It checks the input arguments and sets various tolerances and limits.
%
%  Based on the description for NEINCK.  Part of the modular software system
%  from the appendix of the book "Numerical Methods for Unconstrained
%  Optimization and Nonlinear Equations" by Dennis & Schnabel, 1983.
%
%  Coded in MATLAB by Richard T. Behrens, April 1988.
%

termcode = 0;
dout = din;

% Step 1.
[n,nn] = size(x0);
if ((n < 1) | (nn > 1))
   termcode = -1;
   return
end

% Steps 2 & 3.
[l,m] = size(scale);
if (dout(16) == 2)
   if ((l==n) & (m==2))
      Sx = ones(n,1)./abs(scale(:,1));
      SF = ones(n,1)./abs(scale(:,2));
   else
      dout(16) = 1;
   end
end
if (dout(16) == 0)
   Sx = ones(n,1);
   SF = ones(n,1);
end
if (dout(16) == 1)
   if (any(x0==0))
      x0(find(x0==0)) = ones(sum(x0==0),1);
   end
   if (any(F0==0))
      F0(find(F0==0)) = ones(sum(F0==0),1);
   end
   Sx = ones(n,1)./abs(x0);
   SF = ones(n,1)./abs(F0);
end

% Step 4.
if (dout(12) <= 0)
   dout(13) = eps;
else
   dout(13) = max(eps,10^(-dout(12)));
end
if (dout(13) > .01)
   termcode = -2;
   return
end

% Step 5.
if (dout(2) <= 0)
   dout(2) = 1;    % Default to linesearch.
end
if (((dout(2) == 2) | (dout(2) == 3)) & (dout(7) <= 0))
   dout(7) = -1;
end

% Step 6.
if (dout(6) <= 1)
   dout(6) = 100;  % Default to 100 iteration limit.
end
if (dout(8) <= 0)
   dout(8) = eps ^ (1/3);        % fvectol.
end
if (dout(9) <= 0)
   dout(9) = eps ^ (2/3);        % steptol.
end
if (dout(10) <= 0)
   dout(10) = eps ^ (2/3);       % mintol.
end
if (dout(11) <= 0)
   dout(11) = 1000 * max(norm(Sx .* x0),norm(diag(Sx)));   % maxstep.
end

