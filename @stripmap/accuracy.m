function acc = accuracy(M)
%ACCURACY Apparent accuracy of Schwarz-Christoffel strip map.
%   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel strip
%   map M. The technique used is to compare the differences between
%   successive finite vertices to the integral between the corresponding
%   prevertices, and return the maximum.
%   
%   See also STRIPMAP.

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

% Find ends of strip and renumber to put left end first
endidx(1) = find( isinf(z) & (z < 0) );
endidx(2) = find( isinf(z) & (z > 0) );

% Integrate between consecutive finite pairs along bottom and top
bot = ~imag(z);
bot(endidx) = 0;
top = logical(imag(z));
top(endidx) = 0;
idxbot = find( bot & ~isinf(w) );
idxtop = find( top & ~isinf(w) );

% Two columns hold endpoint indices for integrations
idx = [idxbot(1:end-1) idxbot(2:end)];
idx = [idx ; [idxtop(1:end-1) idxtop(2:end)] ];

% As a final check, integrate once across the strip
[tmp,k] = min(abs( real(z(idxtop)-z(rem(endidx(1),n)+1)) ));
idx = [idx; [rem(endidx(1),n)+1 idxtop(k)]];

I = zeros(size(idx,1),1);
zl = z(idx(:,1));
zr = z(idx(:,2));

% Two-stage integrations (neighboring prevertices)
s2 = (diff(idx,1,2) == 1);
mid = mean([zl(s2) zr(s2)],2);
I(s2) = stquadh(zl(s2),mid,idx(s2,1),z,beta,qdata) ...
    - stquadh(zr(s2),mid,idx(s2,2),z,beta,qdata);

% Three-stage integrations
mid1 = real(zl(~s2)) + i/2;
mid2 = real(zr(~s2)) + i/2;
I(~s2) = stquad(zl(~s2),mid1,idx(~s2,1),z,beta,qdata) ...
    + stquadh(mid1,mid2,zeros(size(mid1)),z,beta,qdata) ...
    - stquad(zr(~s2),mid2,idx(~s2,2),z,beta,qdata);

acc = max(abs( c*I - diff(w(idx).').' ));
