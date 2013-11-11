function [ul,dl] = crspread(u,quadnum,cr,Q)
%CRSPREAD Transform points to every embedding in CR formulation.
%   Each quadrilateral has an associated embedding of the
%   prevertices. These embeddings are linked by Moebius transformations,
%   each of which is well-conditioned.  
%   
%   US = CRSPREAD(U,QUADNUM,CR,Q) assumes that the points of U are given
%   in a single embedding, for quadrilateral QUADNUM. The Moebius
%   transformations are applied recursively so that US(:,K) represents
%   U(K) in all the embeddings. Equivalently, US(QN,:) is the
%   representation of U(:).'  in embedding number QN.
%   
%   [US,DS] = CRSPREAD(U,QUADNUM,CR,Q) also returns the derivatives of
%   the composite transformations to the embeddings.
%       
%   See also CRPARAM, CRGATHER, MOEBIUS.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crspread.m 278 2007-02-28 15:01:27Z driscoll $


n3 = length(cr);

u = u(:).';
ul = zeros(n3,length(u));
ul(quadnum,:) = u;
dl = zeros(n3,length(u));
dl(quadnum,:) = ones(1,length(u));

% Place the quadrilateral prevertices in a rectangle around the origin.
% We will store each rectangle so as to stably find the map between any
% neighboring pair of embeddings.
idx = Q.qlvert(:,quadnum);
r = cr(quadnum);
f1 = sqrt(1/(r+1));
f2 = sqrt(r/(r+1));
zr = zeros(4,n3);
zr(:,quadnum) = f1*[-1;-1;1;1] + i*f2*[-1;1;1;-1];

done = zeros(n3,1);
done(quadnum) = 1;
% Neighbors of quadnum are available
todo = Q.adjacent(:,quadnum);			

while any(~done)
  q = min(find(todo));			% pick an embedding
  r = cr(q);
  % Quadrilateral prevertices for q
  f1 = sqrt(1/(r+1));
  f2 = sqrt(r/(r+1));
  zr(:,q) = f1*[-1;-1;1;1] + i*f2*[-1;1;1;-1];
  % Find a neighbor to quadrilateral q
  qn = min(find(done & Q.adjacent(:,q)));
  % Find the 3 points in common between q and qn
  [i1,i2] = find(Q.qlvert(:,q*ones(4,1))==Q.qlvert(:,qn*ones(4,1))');
  mt = double(moebius(zr(i2,qn),zr(i1,q)));
  ul(q,:) = (mt(2)*ul(qn,:)+mt(1))./(mt(4)*ul(qn,:)+mt(3));
  if nargout > 1
    dl(q,:) = dl(qn,:).*(mt(2)*mt(3)-mt(1)*mt(4))./(mt(4)*ul(qn,:)+mt(3)).^2;
  end
  done(q) = 1;
  todo(q) = 0;
  % Neighbors of q can be done now
  todo = todo | (Q.adjacent(:,q) & ~done);
end
