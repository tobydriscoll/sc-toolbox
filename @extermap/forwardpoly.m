function p = forwardpoly(map)
%   Given an extermap M, FORWARDPOLY(M) returns the polygon that is
%   formed using the prevertices, angles, and quadrature data of that
%   map. If the prevertices were found from the solution of a
%   parameter problem, then the result should agree closely with the
%   original polygon that was supplied.
  
%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: forwardpoly.m 39 1998-07-01 17:40:15Z tad $

z = map.prevertex;
alpha = flipud(angle(polygon(map)));
c = map.constant;

n = length(z);

% Since there is no parameter problem, use high accuracy in quadrature.
qdata = scqdata(1-alpha,16);

% Midpoints of integration
theta = rem(angle(z(n)) + angle(z/z(n))+2*pi,2*pi);
theta(end) = 2*pi;
mid = exp(i*(theta(1:n-1)+theta(2:n))/2);

% Integrations
I = dequad(z(1:n-1),mid,1:n-1,z,1-alpha,qdata) - ...
    dequad(z(2:n),mid,2:n,z,1-alpha,qdata);

% Deduce vertices
w = c*cumsum([0;I]);

p = polygon(w);
