function p = forwardpoly(map)
%   Given a diskmap M, FORWARDPOLY(M) returns the polygon that is
%   formed using the prevertices, angles, and quadrature data of that
%   map. If the prevertices were found from the solution of a
%   parameter problem, then the result should agree closely with the
%   original polygon that was supplied.
  
%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: forwardpoly.m 28 1998-06-22 22:32:30Z tad $

z = map.prevertex;
alpha = angle(polygon(map));
c = map.constant;

n = length(z);

% Since there is no parameter problem, use high accuracy in quadrature.
qdata = sctool.scqdata(alpha-1,16);

w = zeros(n,1);
atinf = (alpha < eps);
w(atinf) = Inf;

% Endpoints of integrations
idx = find(~atinf);
idx = [idx(1:end-1) idx(2:end)];

% Origin is midpoint of every integration
mid = zeros(length(idx),1);

% Integrations
I = dquad(z(idx(:,1)),mid,idx(:,1),z,alpha-1,qdata) - ...
    dquad(z(idx(:,2)),mid,idx(:,2),z,alpha-1,qdata);

% Deduce vertices
w(~atinf) = c*cumsum([0;I]);

p = polygon(w,alpha);
