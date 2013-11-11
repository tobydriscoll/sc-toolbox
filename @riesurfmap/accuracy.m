function acc = accuracy(M)
%ACCURACY Apparent accuracy of a Schwarz-Christoffel map.
%   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel Riemann
%   surface map M. The technique used is to compare the differences between
%   successive finite vertices to the integral between the corresponding
%   prevertices, and return the maximum.
%   
%   See also RIESURFMAP.

%   Copyright 2002 by Toby Driscoll.
%   $Id: accuracy.m 201 2002-09-13 19:46:08Z driscoll $

% If an accuracy has been assigned, don't question it
if ~isempty(M.accuracy)
  acc = M.accuracy;
  return
end

% Get data for low-level functions
p = polygon(M);
w = vertex(p);
beta = angle(p) - 1;
wb = M.branch;
z = M.prevertex;
zb = M.prebranch;
c = M.constant;
qdata = M.qdata;

% Test accuracy by integrating between consecutive finite prevertices, and
% comparing to differences of vertices.

n = length(w);
idx = find(~isinf(w));
wf = w(idx);				% finite vertices

% Two columns hold endpoint indices for integrations.
idx = [idx(1:end) idx([2:end 1])];

% Do the integrations.
I = rsquad(z(idx(:,1)),z(idx(:,2)),idx(:,1),idx(:,2),z,beta,zb,qdata);

acc = max(abs( c*I - diff(wf([1:end 1])) ));
