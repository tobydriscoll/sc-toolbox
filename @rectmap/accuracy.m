function acc = accuracy(M)
%ACCURACY Apparent accuracy of Schwarz-Christoffel rectangle map.
%   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel
%   rectangle map M. The technique used is to compare the differences
%   between successive finite vertices to the integral between the
%   corresponding prevertices, and return the maximum.
%   
%   See also RECTMAP.

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
L = M.stripL;
qdata = M.qdata;

% Renumber to put first corner first
[w,beta,z,corner] = rcorners(w,beta,z);

% Map prevertices to strip
K = max(real(z));
Kp = max(imag(z));
zs = r2strip(z,z,L);
zs = real(zs) + i*round(imag(zs));	% put them *exactly* on edges

% Integrate between consecutive finite pairs on bottom and top
idxbot = find( ~isinf(w(1:corner(3)-1)) );
idxtop = corner(3)-1 + find( ~isinf(w(corner(3):n)) );

% Two columns hold endpoint indices for integrations
idx = [idxbot(1:end-1) idxbot(2:end)];
idx = [idx ; [idxtop(1:end-1) idxtop(2:end)] ];

% Find midpoints. Go into interior of strip to avoid integrating through
% skipped prevertices. 
zz = zs(idx).';
mid = mean(zz).';
mid = real(mid) + i/2;

% As a final check, integrate once across the strip
[tmp,k] = min(abs( real(zs(idxtop)-zs(1)) ));
idx = [idx; [1 idxtop(k)]];
mid(end+1) = mean(zs(idx(end,:)));

% Add in ends of strip
ends = find(diff(imag(zs([1:n 1]))));
zq = [zs(1:ends(1));Inf;zs(ends(1)+1:ends(2));-Inf;zs(ends(2)+1:n)];
bq = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
wq = [w(1:ends(1));NaN;w(ends(1)+1:ends(2));NaN;w(ends(2)+1:n)];
% Extend qdat with useless columns at ends
j = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1];
qdata = qdata(:,[j j+n+1]);

% Do the integrations
zleft = zs(idx(:,1));
zright = zs(idx(:,2));
idx = idx + (idx > ends(1)) + (idx > ends(2));
I = stquad(zleft,mid,idx(:,1),zq,bq,qdata) - ...
    stquad(zright,mid,idx(:,2),zq,bq,qdata);

acc = max(abs( c*I - diff(wq(idx).').' ));
