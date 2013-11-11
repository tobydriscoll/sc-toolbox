function [hits,loc] = intersect(p,endpt,tol)
%INTERSECT  Find intesection of segments with polygon sides.
%
%   S = INTERSECT(P,ENDPT) checks for intesections between sides of the
%   polygon P and the line segments whose endpoints are given in the
%   complex M by 2 matrix ENDPT. If P has N sides, on return S is an
%   M by N logical matrix with nonzeros at locations indicating
%   intersection. 
%
%   INTERSECT(P,ENDPT,TOL) requires that the intersection take place
%   more than TOL away (relatively) from the segments' endpoints. By
%   default TOL=EPS. To test truly closed segments, use
%   INTERSECT(P,ENDPT,0); however, this is a poorly conditioned
%   problem.

%   Copyright 2002 by Toby Driscoll. 
%   $Id: intersect.m 189 2002-09-10 17:18:17Z driscoll $

n = length(p);
if nargin < 3
  tol = eps;
end

m = size(endpt,1);
w = vertex(p);
beta = angle(p)-1;

% Where are the slits?
isslit = abs(beta-1) < 2*eps;
isslit = isslit | isslit([2:n 1]);

% Find two consecutive finite vertices.
dw = diff( w([1:n 1]) );
K = min( find( ~isinf(dw) ) );
% Arguments of polygon sides.
argw = ones(n,1);
argw([K:n 1:K-1]) = cumsum( [angle(dw(K));-pi*beta([K+1:n 1:K-1])] );

% Check each side. Solve for two parameters and check their ranges.
hits = logical(zeros(m,n));
loc = repmat(NaN,[m n]);
for k = 1:n
  tangent = exp(i*argw(k));
  if ~isinf(w(k))
    wk = w(k);
    s1max = abs( w(rem(k,n)+1)-w(k) );  % parameter in [0,s1max]
  else
    % Start from next vertex and work back.
    wk = w(rem(k,n)+1);
    tangent = -tangent;
    s1max = Inf;
  end
  A(:,1) = [ real(tangent); imag(tangent) ];

  % Loop over the segments to be tested. The alternative is to solve a
  % block 2x2 diagonal matrix, but any collinear cases would ruin the
  % whole batch.
  for j = 1:m
    e1e2 = endpt(j,2) - endpt(j,1);
    A(:,2) = -[ real(e1e2); imag(e1e2) ];
    if rcond(A) < 2*eps
      % Segments are parallel. Check for collinearity using rotation.
      e2 = (endpt(j,2)-wk) / tangent;
      e1 = (endpt(j,1)-wk) / tangent;
      if abs(imag(e1)) < 2*eps
        % Check for overlapping.
        x1 = min( real([e1 e2]) );
        x2 = max( real([e1 e2]) );
        % Do these values straddle either of the side's endpoints?
        if (x2 >= tol) & (x1 <= s1max-tol)
          hits(j,k) = 1;
          loc(j,k) = wk;  % pick a place
        end
      end
    else
      % Generic case. Find intersection.
      delta = endpt(j,1) - wk;
      s = A \ [real(delta);imag(delta)];
      % Check parameter ranges.
      if s(1)>=-eps & s(1)<=s1max+eps & s(2)>=tol & s(2)<=1-tol
        % If an end of the segment lies on a slit side, check for
        % interior vs. exterior.
        if isslit(k) & (abs(s(2)) < 10*eps)
          normal = i*tangent; 
          if real( conj(e1e2)*normal ) < 0, break, end
        elseif isslit(k) & (abs(s(2)-1) < 10*eps)
          normal = i*tangent;
          if real( conj(e1e2)*normal ) > 0, break, end
        end
        hits(j,k) = 1;
        loc(j,k) = wk + s(1)*tangent;
      end
    end
  end
end
