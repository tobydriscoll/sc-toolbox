function p = forwardpoly(map)
%   Given an hplmap M, FORWARDPOLY(M) returns the polygon that is
%   formed using the prevertices, angles, and quadrature data of that
%   map. If the prevertices were found from the solution of a
%   parameter problem, then the result should agree closely with the
%   original polygon that was supplied.
  
%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: forwardpoly.m 37 1998-06-30 00:22:40Z tad $

z = map.prevertex;
alpha = angle(polygon(map));
c = map.constant;

n = length(z);

% Since there is no parameter problem, use high accuracy in quadrature.
qdata = scqdata(alpha(1:n-1)-1,16);

w = zeros(n,1);
atinf = (alpha < eps);
w(atinf) = Inf;

% Endpoints of integrations. Because the last prevertex is at Inf, we
% shouldn't try to integrate there.
idx = find(~atinf);
if idx(end)==n, idx(end) = []; end
endpt = [idx(1:end-1) idx(2:end)];

% Midpoints are in upper half-plane. Always make 45 degrees with real line.
mid = mean(z(endpt),2) + 1i*diff(z(endpt),1,2)/2;

% Integrations
I = hpquad(z(endpt(:,1)),mid,endpt(:,1),z(1:n-1),alpha(1:n-1)-1,qdata) - ...
    hpquad(z(endpt(:,2)),mid,endpt(:,2),z(1:n-1),alpha(1:n-1)-1,qdata);

% Deduce vertices
w(idx) = c*cumsum([0;I]);

% Get the last vertex via intersection
if alpha(n) > 0
  if abs(alpha(n)-1) < 5*eps || abs(alpha(n)-2) < 5*eps
    error(['Cannot deduce last vertex when its adjacent sides are' ...
          ' collinear.'])
  elseif any(atinf([1 2 n-1]))
    error('Vertices 1, 2, and end-1 must be finite.')
  else
    % Here's the direction from w(1)
    d1 = (w(2)-w(1))*exp(1i*pi*alpha(1));
    % Get the direction from w(n-1)
    d2 = angle(w(2)-w(1)) + sum(pi*(1-alpha(2:n-1)));
    d2 = exp(1i*d2);
    b = w(n-1) - w(1);
    s = [real([d1 -d2]);imag([d1 -d2])]\[real(b);imag(b)];
    w(n) = w(1) + s(1)*d1;
  end
end
    
p = polygon(w,alpha);
