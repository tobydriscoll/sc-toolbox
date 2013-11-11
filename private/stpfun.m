function F = stpfun(y,fdat)
%   Returns residual for solution of nonlinear equations. 
%   $Id: stpfun.m 62 1999-01-29 00:56:34Z tad $

[n,nb,beta,nmlen,left,right,cmplx,qdat] = deal(fdat{:});

% In this function, n refers to the number of FINITE prevertices.

% Transform y (unconstr. vars) to z (actual params)
z = zeros(n,1);
z(2:nb) = cumsum(exp(y(1:nb-1)));
z(nb+1:n) = i+cumsum([y(nb);-exp(y(nb+1:n-1))]);

% Compute the integrals 
zleft = z(left);
zright = z(right);
mid = mean([zleft.' ; zright.']).';
c2 = cmplx;
c2(2) = 0;
mid(c2) = mid(c2) - sign(left(c2)-nb)*i/2;

% Add ends of strip to z, and modify singularity indices
zs = [-Inf;z(1:nb);Inf;z(nb+1:n)];
left = left + 1 + (left > nb);
right = right + 1 + (right > nb);

% Do those staying on a side
ints = zeros(n-1,1);
c2(1) = 1;
id = ~c2;
ints(id) = stquadh(zleft(id),mid(id),left(id),zs,beta,qdat) - ...
    stquadh(zright(id),mid(id),right(id),zs,beta,qdat);

% For the rest, go to the strip middle, across, and back to the side
z1 = real(zleft(c2)) + i/2;
z2 = real(zright(c2)) + i/2;
id = ~id;
ints(id) = stquad(zleft(id),z1,left(id),zs,beta,qdat);
ints(id) = ints(id) + stquadh(z1,z2,zeros(size(z1)),zs,beta,qdat);
ints(id) = ints(id) - stquad(zright(id),z2,right(id),zs,beta,qdat);

absval = abs(ints(~cmplx)); 		% absval(1) = abs(ints(1))
if ~absval(1)
  rat1 = 0;
  rat2 = 0;
else
  rat1 = absval(2:length(absval))/absval(1);
  rat2 = ints(cmplx)/ints(1);
end

if any([rat1;rat2]==0) | any(isnan([rat1;rat2])) | any(isinf([rat1;rat2]))
  % Singularities were too crowded. 
  warning('Severe crowding')
end

% Compute nonlinear equation residual values.
cmplx2 = cmplx(2:length(cmplx));
if ~isempty(rat1)
    F1 = log( rat1 ./ nmlen(~cmplx2) );
else
  F1 = [];
end
if ~isempty(rat2)
  F2 = log( rat2 ./ nmlen(cmplx2) );
else
  F2 = [];
end
F = [F1;real(F2);imag(F2)];

