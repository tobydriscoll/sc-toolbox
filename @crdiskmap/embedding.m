function [z,c,a] = embedding(M,j)
%EMBEDDING Embeddings of prevertices.
%   The CR formulation does not use a single embedding (configuration)
%   of the prevertices. Rather, for an n-gon it uses n-3 of them,
%   related by Mobius transformations. Each is accurate when evaluating
%   the S--C integral for points inside a quadrilateral in the domain
%   decomposition.
%   
%   Z = EMBEDDING(M,J) returns the prevertices for the embedding(s)
%   numbered J. If J has more than one entry, Z will be a cell
%   array. EMBEDDING(M) returns all embeddings.
%   
%   [Z,C,A] = EMBEDDING(M,J) also returns the multiplicative and
%   additive constant(s) needed to properly map to POLYGON(M) when
%   evaluating the S--C integral from the origin.
% 
%   DISKMAP(Z,C) will create an ordinary S--C disk map for an embedding.
%   
%   See also CRDISKMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: embedding.m 7 1998-05-10 04:37:19Z tad $

% One arg: Do them all
if nargin == 1
  j = 1:length(polygon(M))-3;
end

% Empty return
if isempty(j)
  z = [];
  c = [];
  a = [];
  return
end

% Pre-allocate
z = cell(size(j));
c = cell(size(j));
a = cell(size(j));

% Fill them
for k = 1:length(j(:))
  z{k} = crembed(M.crossratio,M.qlgraph,j(k));
  c{k} = M.affine(j(k),1);
  a{k} = M.affine(j(k),2);
end

% Don't return cells for just one
if isscalar(j(:))
  z = z{1};
  c = c{1};
  a = a{1};
end
