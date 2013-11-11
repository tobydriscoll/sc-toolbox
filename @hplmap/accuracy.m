function acc = accuracy(M)
%ACCURACY Apparent accuracy of Schwarz-Christoffel half-plane map.
%   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel
%   half-plane map M. The technique used is to compare the differences
%   between successive finite vertices to the integral between the
%   corresponding prevertices, and return the maximum.
%   
%   See also HPLMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: accuracy.m 7 1998-05-10 04:37:19Z tad $

% If an accuracy has been assigned, don't question it
if ~isempty(M.accuracy)
  acc = M.accuracy;
  return
end

% Get data for low-level functions
p = polygon(M);
w = vertex(p);
n = length(w);
beta = angle(p) - 1;
z = M.prevertex;
c = M.constant;
qdata = M.qdata;

% Test accuracy by integrating between consecutive finite prevertices, and
% comparing to differences of vertices.
n = length(w);	
idx = find(~isinf(w(1:n-1)));		% exclude last prevert, at Inf

% Two columns hold endpoint indices for integrations
idx = [idx(1:end-1) idx(2:end)];

% Find midpoints. Go into upper half-plane to avoid integrating through
% skipped prevertices.
zz = z(idx).';
mid = mean(zz).';
mid = mid + i*abs(diff(zz).')/2;

% Do the integrations
I = hpquad(z(idx(:,1)),mid,idx(:,1),z(1:n-1),beta(1:n-1),qdata) - ...
    hpquad(z(idx(:,2)),mid,idx(:,2),z(1:n-1),beta(1:n-1),qdata);

acc = max(abs( c*I - diff(w(idx).').' ));
