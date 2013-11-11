function F = hppfun(y,fdat)
%   Returns residual for solution of nonlinear equations. 
%   $Id: hppfun.m 60 1999-01-29 00:49:09Z tad $

[n,beta,nmlen,left,right,cmplx,qdat] = deal(fdat{:});

% Transform y (unconstr. vars) to z (prevertices)
cp = cumprod([1;exp(-y)]);
z = [0;cumsum(cp)] - [flipud(cumsum(flipud(cp)));0];
z = z/z(n-1);

% Compute the integrals 
zleft = z(left);
zright = z(right);
mid = mean([zleft.' ; zright.']).';
% For integrals between non-adjacent singularities, choose intermediate
% points in the upper half-plane.
mid(cmplx) = mid(cmplx) + i*(zright(cmplx)-zleft(cmplx))/2;
ints = hpquad(zleft,mid,left,z,beta,qdat) - ...
    hpquad(zright,mid,right,z,beta,qdat);

if any(ints==0)
  % Singularities were too crowded in practice.
  warning('Severe crowding')
end

% Compute nonlinear equation residual values.
F1 = abs(ints(~cmplx));		% F1(1) = abs(ints(1))
F1 = F1(2:end)/F1(1);
F2 = ints(cmplx)/ints(1);
F = [F1;real(F2);imag(F2)] - nmlen;

%%ints = ints/ints(1);
%%w = [0;cumsum([1;ints(2:end)])];
%%set(findobj(0,'tag','polydata'),'xd',real(w),'yd',imag(w))
%%drawnow
