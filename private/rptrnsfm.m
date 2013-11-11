function z = rptrnsfm(y,cnr)
%RPTRNSFM (not intended for calling directly by the user)
%   Transform optimization vars to prevertices for rectangle parameter
%   problem.

%       Copyright 1997 by Toby Driscoll. Last updated 05/06/97.

n = length(y)+3;
z = zeros(n,1);

% Fill interior of "long edges" first
z(cnr(1)+1:cnr(2)-1) = cumsum(exp(y(cnr(1):cnr(2)-2)));
z(cnr(4)-1:-1:cnr(3)+1) = i + cumsum(exp(y(cnr(4)-3:-1:cnr(3)-1)));

% Find L
xr = real( z([cnr(2)-1,cnr(3)+1]) );
z(cnr(2)) = mean(xr)+sqrt(diff(xr/2)^2+exp(2*y(cnr(2)-1)));
z(cnr(3)) = i + z(cnr(2));
z(cnr(4)) = i;

% Now, fill in "short edges"
cp = cumprod([1;exp(-y(cnr(2):cnr(3)-2))]);
x = [0;cumsum(cp)] - [flipud(cumsum(flipud(cp)));0];
x = x(2:end-1)/x(end);
mask = abs(x) < eps;
u = x;
u(~mask) = log( x(~mask) ) / pi;
u(mask) = -z(cnr(2))/eps;
z(cnr(2)+1:cnr(3)-1) = i*imag(u) + real(z(cnr(2))) - real(u);

idx = [cnr(4)-2:n-3 1:cnr(1)-1];
cp = cumprod([1;exp(-y(idx))]);
x = [0;cumsum(cp)] - [flipud(cumsum(flipud(cp)));0];
x = x(2:end-1)/x(end);
mask = abs(x) < eps;
u = x;
u(~mask) = log( x(~mask) ) / pi;
u(mask) = -z(cnr(2))/eps;
z([cnr(4)+1:n 1:cnr(1)-1]) = u;
