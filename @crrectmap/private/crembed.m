function z = crembed(cr,Q,qnum)
%CREMBED Embed prevertices for given crossratios.
%   CREMBED(CR,Q,QNUM) embeds prevertices of given crossratios CR and
%   quadrilateral graph information Q (see QLGRAPH). The embedding is
%   constructed so that the prevertices of quadrilateral QNUM are in a
%   rectangle and do not crowd up against each other or any other
%   prevertices (as long as no crossratios are very far from unity). The
%   resulting S-C map should be accurate for that rectangle, inculding
%   the four prevertices of quadrilateral QNUM.
    
%   Copyright 1998 by Toby Driscoll.
%   $Id: crembed.m 7 1998-05-10 04:37:19Z tad $

%   This function is adapted from a C routine written by Stephen Vavasis.

n = length(cr)+3;
z = zeros(n,1);

% Place the quadrilateral spanned by the initial diagonal
% at a rectangle around the origin.
r = cr(qnum);
f1 = sqrt(1/(r+1));
f2 = sqrt(r/(r+1));
idx = Q.qlvert(:,qnum);			% quadrilateral vertices
z(idx) = f1*[-1;-1;1;1] + 1i*f2*[-1;1;1;-1];

% Set up "already visited" lists
vtxdone = zeros(n,1);
vtxdone(idx) = ones(4,1);
edgedone = zeros(2*n-3,1);
edgedone(qnum) = 1;
% Mark boundary edges as "done" so they will be ignored
edgedone(n-2:2*n-3) = ones(n,1);

% Set up "ready to do" list. There does not need to be a strict ordering as
% in a queue or stack; any ready edge can go next.
edgetodo = zeros(2*n-3,1);
% The edges of the initial quadrilateral are ready to go
edgetodo(Q.qledge(:,qnum)) = ~edgedone(Q.qledge(:,qnum));

% Begin iteration
while any(edgetodo)
  e = min(find(edgetodo));
  idx = Q.qlvert(:,e);
  % If necessary, renumber so that z(idx(4)) is to be determined
  if ~vtxdone(idx(2))
    idx = idx([3:4,1:2]);
  end
  
  % Work around divide by zero. This will have no dire effects,
  % since there is no guarantee that the prevertices not in qnum will be
  % separated.
  if abs(z(idx(2))-z(idx(1))) < 5*eps
    z(idx(4)) = z(idx(2));
  else
    % Place z(idx(4))
    fac = -cr(e)*(z(idx(3)) - z(idx(2))) / (z(idx(2)) - z(idx(1)));
    z(idx(4)) = (z(idx(3)) + fac*z(idx(1))) / (1 + fac);
    % Keep the vertices exactly on the unit circle (fix roundoff)
    z(idx(4)) = sign(z(idx(4)));
  end
  
  vtxdone(idx(4)) = 1;
  edgedone(e) = 1;
  edgetodo(e) = 0;
  % Ready to do neighboring edges that still have not been visited.
  edgetodo(Q.qledge(:,e)) = ~edgedone(Q.qledge(:,e)); 
    
end    
