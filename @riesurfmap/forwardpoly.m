function p = forwardpoly(map)
%   Given a riesurfmap F, FORWARDPOLY(F) returns the polygon that is
%   formed using the prevertices, angles, and quadrature data of that
%   map. If the prevertices were found from the solution of a
%   parameter problem, then the result should agree closely with the
%   original polygon that was supplied.
  
%   Copyright 2002 by Toby Driscoll.
%   $Id: forwardpoly.m 203 2002-09-13 20:09:37Z driscoll $

z = map.prevertex;
alpha = angle(polygon(map));
c = map.constant;
zb = map.prebranch;

n = length(z);

% Since there is no parameter problem, use high accuracy in quadrature.
qdata = scqdata(alpha-1,16);

w = zeros(n,1);
atinf = (alpha < eps);
w(atinf) = Inf;

% Endpoints of integrations
idx = find(~atinf);
idx = [idx(1:end-1) idx(2:end)];

% Integrations
I = rsquad(z(idx(:,1)),z(idx(:,2)),idx(:,1),idx(:,2),z,alpha-1,zb,qdata);

% Deduce vertices
w(~atinf) = c*cumsum([0;I]);

p = polygon(w,alpha);
