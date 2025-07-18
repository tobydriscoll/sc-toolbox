function acc = accuracy(M)
%ACCURACY Apparent accuracy of Schwarz-Christoffel disk map.
%   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel disk
%   map M. The technique used is to compare the differences between
%   successive finite vertices to the integral between the corresponding
%   prevertices, and return the maximum.
%   
%   See also DISKMAP.

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
beta = angle(p) - 1;
z = M.prevertex;
c = M.constant;
qdata = M.qdata;

% Test accuracy by integrating between consecutive finite prevertices, and
% comparing to differences of vertices.

idx = find(~isinf(w));
wf = w(idx);				% finite vertices

% Two columns hold endpoint indices for integrations
idx = [idx(1:end) idx([2:end 1])];

% Always use center as the integration midpoint
%dtheta = mod(angle(z(idx(:,2))./z(idx(:,1))),2*pi);
%mid = z(idx(:,1)).*exp(1i*dtheta/2);
mid = zeros(length(idx),1);

% Do the integrations
I = dquad(z(idx(:,1)),mid,idx(:,1),z,beta,qdata) - ...
    dquad(z(idx(:,2)),mid,idx(:,2),z,beta,qdata);

acc = max(abs( c*I - diff(wf([1:end 1])) ));
