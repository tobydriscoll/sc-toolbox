function [edge,triedge,edgetri] = crcdt(w,edge,triedge,edgetri)
%CRCDT  Constrained Delaunay triangulation of a polygon.
%   [E,TE,ET] = CRCDT(W,E,TE,ET) computes a constrained Delaunay
%   triangulation of a polygon given any initial triangulation. The
%   parameters E, TE, and ET describe the triangulation as in POLYTRI.
       
%   Copyright 1998 by Toby Driscoll.
%   $Id: crcdt.m 7 1998-05-10 04:37:19Z tad $

numedge = size(edge,2);
interior = (edgetri(2,:)~=0);
% done marks those that are known to be correct. 
done = ones(numedge,1);	
done(interior) = zeros(sum(interior),1); % boundaries are fixed
quadvtx = zeros(4,1);

while any(~done)
  e = find(~done, 1 );

  % Get the 2 triangles in which edge e participates
  t1 = triedge(:,edgetri(1,e));
  t2 = triedge(:,edgetri(2,e));

  % Find the quadrilateral of which e is the diagonal
  quadvtx([1,3]) = edge(:,e);
  e1 = edge(:,t1); e1 = sort(e1(:)); e1 = e1(1:2:5);
  e2 = edge(:,t2); e2 = sort(e2(:)); e2 = e2(1:2:5);
  quadvtx(2) = e1((e1~=quadvtx(1))&(e1~=quadvtx(3)));
  quadvtx(4) = e2((e2~=quadvtx(1))&(e2~=quadvtx(3)));

  % Find the angles where diagonal meets quadrilateral
  edgew = diff(w(edge(:,e)));
  alpha(1) = angle(edgew/diff(w(quadvtx(1:2))));
  alpha(2) = angle(diff(w(quadvtx([1,4])))/edgew);
  alpha(3) = angle(-edgew/diff(w(quadvtx(3:4))));
  alpha(4) = angle(-diff(w(quadvtx([3,2])))/edgew);
  alpha = abs(alpha);
  
  % Flip if sum of angles is < pi
  if sum(alpha) < pi
    
    % Change endpts of edge e
    edge(:,e) = quadvtx([2,4]);
    % Find a quadrilateral side in each triangle
    i1 = rem(find(t1==e),3)+1;
    i2 = rem(find(t2==e),3)+1;
    % They must be opposites in the quadrilateral (no common endpts)
    if any(any( (edge(:,t1([i1 i1])) - edge(:,t2([i2 i2]))')==0 ))
      i2 = rem(find(t2==e)+1,3)+1;
    end
    % Change triangle/edge assignments
    triedge(i1,edgetri(1,e)) = t2(i2);
    edgetri(edgetri(:,t2(i2))==edgetri(2,e),t2(i2)) = edgetri(1,e);
    triedge(i2,edgetri(2,e)) = t1(i1);
    edgetri(edgetri(:,t1(i1))==edgetri(1,e),t1(i1)) = edgetri(2,e);
    
    % The (non-bdy) edges of the triangles must be reconsidered
    done(t1) = ~interior(t1);
    done(t2) = ~interior(t2);

  end
  % Edge e is done, unless a neighbor resets it
  done(e) = 1;
      
end
