function [edge,triedge,edgetri] = crtriang(w)
%CRTRIANG Triangulate a polygon.
%   [E,TE,ET] = CRTRIANG(W) triangulates a polygon at the N vertices
%   W. On exit, E is a 2 x (2N-3) matrix of edge endpoint indices
%   (interior edges in columns 1:N-3), TE is a 3 x (N-2) matrix of
%   triangle edge indices, and ET is a 2 x (2N-3) matrix of triangle
%   membership indices for the edges. If edge k is a member of just one
%   triangle (a boundary edge), then ET(2,k)=0.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crtriang.m 7 1998-05-10 04:37:19Z tad $

import sctool.*

% Uses a stack-based approach. The recursive version is simpler but prone to
% causing memory errors.

W = w;
N = length(W);

% Initialize outputs.
edge = zeros(2,2*N-3);
triedge = zeros(3,N-2);
edgetri = zeros(2,2*N-3);

% Enter original boundary edges at the end of edge.
edge(1,N-3+(1:N)) = 1:N;
edge(2,N-3+(1:N)) = [2:N 1];

% Each row is a (potential) stack entry, consisting of N entries of 0-1
% indices. These indices describe a subdivision of the original polygon. We
% preallocate memory to avoid repeated resizing that could bog down. To
% avoid hogging too much memory if N is large, we will hope that the number
% of polygon subdivisions does not exceed 300.
maxstack = min(N,300);
stack = NaN*zeros(maxstack,N);	

% Intialize stack
stack(1,:) = ones(1,N);
stackptr = 1;

while stackptr > 0
  idx = find(stack(stackptr,:));
  stackptr = stackptr - 1;

  w = W(idx);
  n = length(w);
  
  if n==3	
    % Base case: polygon is a triangle

    % Assign the triangle a new number
    tnum = min(find(any(triedge==0)));

    e = zeros(2,3);
    e(:) = idx([1 2 2 3 3 1]);		% edge vertices
    e = sort(e);
    de = diff(e);
    for j=1:3
      % Find reference # of triangle edge j
      if de(j)==1
	% Adjacent vertices
	enum = N-3 + e(1,j);
      elseif de(j)==N-1
	% Vertices are [1 N]
	enum = 2*N-3;
      else
	% Find the edge # among the first N-3
	enum = find(all( edge(:,1:N-3) - e(:,j)*ones(1,N-3) == 0 ));
      end
      % Assign edge # to triangle and triangle # to edge
      triedge(j,tnum) = enum;
      edgetri(min(find(edgetri(:,enum)==0)),enum) = tnum;
    end

  else
    
    % More than 3 vertices, so add an edge and subdivide

    % Start with sharpest outward point in polygon
    beta = scangle(w);
    [junk,j] = min(beta);
  
    % Trial triangle has vertices at j and its neighbors
    jm1 = rem(j+n-2,n)+1;
    jp1 = rem(j,n)+1;
    t = logical(zeros(n,1));
    t([jm1;j;jp1]) = ones(3,1);

    % Determine which remaining vertices lie in trial triangle
    triangle = w(t);
    inside = zeros(n,1);
    inside(~t) = abs(isinpoly(w(~t),triangle));
    % Borderline cases: on a triangle edge.
    for k = 1:3
      d = ones(1,n);
      d(~t) = crpsdist(triangle(rem(k+(0:1),3)+1),w(~t));
      for p = find(d < eps)
	inside(p) = 0;
	% If vertex is coincident with a triangle vertex, not "inside"
	if min(abs(w(p)-triangle)) > eps
	  % Polygon slits/cracks need special care
	  % If inward perturbation goes into triangle, it's "inside"
	  s = exp(i*(angle(w(rem(p-2+n,n)+1)-w(p)) - pi*(beta(p)+1)/2));
	  if isinpoly(w(p)+1e-13*abs(w(p))*s,triangle)
	    inside(p) = 1;
	  end
	end
      end
    end  
  
    inside = find(inside);
    if isempty(inside)
      % The trial triangle is OK; use its new edge 
      e = sort([jm1;jp1]);
    else

      % Edge must be drawn from w(j) to a vertex inside t
      
      if length(inside) > 1
	% thanks to S. Vavasis for following
	% Vertices must be "visible" to each other
	% Find angle between forward side and connecting ray
	fwd = rem(inside,n) + 1;
	ang1 = angle( (w(j)-w(inside)) ./ (w(fwd)-w(inside)) );
	ang1 = rem(ang1/pi + 2, 2);
	ang2 = angle( (w(inside)-w(j)) ./ (w(jp1)-w(j)) );
	ang2 = rem(ang2/pi + 2, 2);
	% Detect visibility by requiring ang to be less than interior angle
	vis = find( (ang1 < beta(inside)+1) & (ang2 < beta(j)+1) );
	
	% Find a line through w(j) outside the polygon
	dw = [w(jp1)-w(j) w(j)-w(jm1)];
	theta = angle(dw(1)) + angle(dw(2)/dw(1))/2;
	
	% Find nearest visible point to that line
	wvis = w(inside(vis)) - w(j);
	D = abs(wvis - real(wvis*exp(-i*theta))*exp(i*theta));
	[junk,k] = min(D);
      else
	vis = 1;
	k = 1;
      end
      
      e = sort([inside(vis(k)); j]);
      
    end
    % Assign the next available edge #
    enum = min(find(any(edge==0)));
    edge(:,enum) = idx(e)';
    
    % Indices of subdivided pieces
    i1 = e(1):e(2);
    i2 = [e(2):n 1:e(1)];

    % If stack will overflow, allocate a big chunk of memory
    if stackptr > maxstack-2
      addlen = min(maxstack,N-maxstack);
      stack = [stack;zeros(addlen,N)];
      maxstack = maxstack + addlen;
    end
      
    % Put the two new pieces on the stack
    stackptr = stackptr + 1;
    stack(stackptr,:) = zeros(1,N);
    stack(stackptr,idx(i1)) = ones(1,length(i1));
    stackptr = stackptr + 1;
    stack(stackptr,:) = zeros(1,N);
    stack(stackptr,idx(i2)) = ones(1,length(i2));
    
  end

end
