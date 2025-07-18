function Q = crqgraph(w,edge,triedge,edgetri)
%CRQGRAPH Quadrilateral graph of a triangulation.
%   Q = CRQGRAPH(W,E,TE,ET) constructs the "quadrilateral graph" of a
%   polygon triangulation. TE, and ET are as in CRTRIANG.
%   
%   On return, Q is a structure. Q.edge is 2x(2n-3), each column being
%   the endpoint indices of an edge of the triangulation. The n-3
%   interior edges (diagonals) come first. Q.qledge is 4x(n-3), column k
%   being the indices of the four edges (in counterclockwise order) of
%   the quadrilateral of which edge k is a diagonal. Q.qlvert is
%   4x(n-3), with column k being the indices of the four vertices (in
%   clockwise order) of that quadrilateral. An endpoint of the diagonal
%   of the quadrilateral will be the first entry of the column.
%   Q.adjacent is an (n-3)x(n-3) logical matrix; it indicates which
%   quadrilaterals share a common triangle.
%   
%   If called with only one input argument, CRQGRAPH will call CRTRIANG
%   and CRCDT to create a Delaunay triangulation.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crqgraph.m 7 1998-05-10 04:37:19Z tad $

import sctool.*
% Find triangulation if requested
if nargin==1
  [edge,triedge,edgetri] = crtriang(w);
  [edge,triedge,edgetri] = crcdt(w,edge,triedge,edgetri);
end

% Get # of quadrilaterals
N = size(edge,2);
n3 = (N+3)/2 - 3;

% Now construct the graph
qlvert = zeros(4,n3);
qledge = zeros(4,n3);
T = zeros(n3);
for e=1:n3
  % Find the triangles that include edge e
  t1 = triedge(:,edgetri(1,e));
  t2 = triedge(:,edgetri(2,e));

  % Re-order triangles so that edge e is first
  i1 = find(t1==e);
  t1 = t1([i1:3,1:i1-1]);  
  i2 = find(t2==e);
  t2 = t2([i2:3,1:i2-1]);
  
  % Ensure ccw ordering of edges
  we = zeros(2,3);
  we(:) = w(edge(:,t1));
  if sum(scangle(mean(we))) > 0
    t1(2:3) = t1([3 2]);
  end
  we(:) = w(edge(:,t2));
  if sum(scangle(mean(we))) > 0
    t2(2:3) = t2([3 2]);
  end
    
  % Read off quadrilateral edges
  qledge(:,e) = [t1(2:3);t2(2:3)];
  
  % Any ql edge that is also a diagonal, is adjacent
  t = qledge(qledge(:,e)<=n3,e);
  T(t,e) = ones(length(t),1);
  T(e,t) = ones(1,length(t));
  
  % Extract vertices
  % An endpoint of e must be first (and third)
  qlvert([1,3],e) = edge(:,e);
  % Find unique vertices of t1 and t2 to get other ql vertices
  e1 = edge(:,t1); 
  e1 = sort(e1(:)); 
  e1 = e1(1:2:5);
  e2 = edge(:,t2); 
  e2 = sort(e2(:)); 
  e2 = e2(1:2:5);
  qlvert(2,e) = e1((e1~=qlvert(1,e))&(e1~=qlvert(3,e)));
  qlvert(4,e) = e2((e2~=qlvert(1,e))&(e2~=qlvert(3,e)));
  % Reverse ordering if necessary to clockwise
  if abs(sum(scangle(w(qlvert(:,e))))-2) > 1e-8
    qlvert([2,4],e) = qlvert([4,2],e);
  end
  
end

% Output form
Q.edge = edge;
Q.qledge = qledge;
Q.qlvert = qlvert;
Q.adjacent = logical(T);
