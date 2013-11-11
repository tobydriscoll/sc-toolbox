function wc = center(M,wc)
%CENTER Conformal center of Schwarz-Christoffel disk map.
%   CENTER(M) returns the conformal center (image of 0) of the
%   Schwarz-Christoffel crossratio disk map represented by M.
%   
%   CENTER(M,WC) computes a map conformally equivalent to M but with
%   conformal center WC (provided WC is inside the polygon of M), and
%   returns the new map. If WC is empty, you will be asked to select it
%   graphically. 
%   
%   See also CRDISKMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: center.m 7 1998-05-10 04:37:19Z tad $

if nargin == 1
  wc = M.center{1};
else
  p = polygon(M);
  cr = M.crossratio;
  w = vertex(p);
  beta = angle(p) - 1;
  
  if isempty(wc)
    [wcfix,wc] = crfixwc(w,beta,cr,M.affine,M.qlgraph);
  else
    wcfix = crfixwc(w,beta,cr,M.affine,M.qlgraph,wc);
  end
  
  M.center = {wc,wcfix};
    
  wc = M;
end
