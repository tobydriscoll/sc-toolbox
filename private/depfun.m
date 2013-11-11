function F = depfun(y,fdat)
%   Returns residual for solution of nonlinear equations.
%   $Id: depfun.m 268 2003-05-08 18:08:25Z driscoll $

[n,beta,nmlen,qdat] = deal(fdat{:});

% Transform y (unconstr. vars) to z (prevertices)
cs = cumsum(cumprod([1;exp(-y)]));
theta = 2*pi*cs(1:n-1)/cs(n);
z = ones(n,1);
z(1:n-1) = exp(i*theta);

% Compute the integrals 
mid = exp(i*(theta(1:n-2)+theta(2:n-1))/2);

% We can use the same quadrature as for the interior map, because the abs
% value of the integrand on the unit circle is not affected by the z^{-2}
% term. 
ints = dabsquad(z(1:n-2),mid,1:n-2,z,beta,qdat) + ...
    dabsquad(z(2:n-1),mid,2:n-1,z,beta,qdat);

if any(ints==0)
  % Singularities were too crowded in practice.
  warning('Severe crowding')
end

% Compute equation residual values.
if n > 3
  F = abs(ints(2:n-2))/abs(ints(1)) - nmlen;  
else
  F = [];
end

% Compute residue.
res = -sum(beta./z)/ints(1);

F = [F;real(res);imag(res)];

