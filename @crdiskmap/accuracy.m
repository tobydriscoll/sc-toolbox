function acc = accuracy(M)
%ACCURACY Apparent accuracy of Schwarz-Christoffel cross-ratio disk map.
%   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel CR disk
%   map M. The technique used is to compare the cross-ratios of the actual
%   polygon image with those of the target polygon, and return the maximum.
%   
%   See also CRDISKMAP.

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
cr = M.crossratio;
aff = M.affine;
Q = M.qlgraph;
qdata = M.qdata;

n = length(w);

% Crossratios of target polygon
crtarget = crossrat(w,Q);

% Actual crossratios
crimage = zeros(n-3,1);			% image vertex crossratios

% Compute crossratio for each image quadrilateral
for k = 1:n-3
  prever = crembed(cr,Q,k);
  wq = -crquad(prever(Q.qlvert(:,k)),Q.qlvert(:,k),prever,beta,qdata);
  crimage(k) = (wq(2)-wq(1))*(wq(4)-wq(3))/((wq(3)-wq(2))*(wq(1)-wq(4)));
end

% Compare them
acc = max(abs(crimage-crtarget));
