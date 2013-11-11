function wc = center(map,wc)
%CENTER Conformal center of Schwarz-Christoffel disk map.
%   CENTER(M) returns the conformal center (image of 0) of the
%   Schwarz-Christoffel disk map represented by M.
%   
%   CENTER(M,WC) computes a map conformally equivalent to M but with
%   conformal center WC (provided WC is inside the polygon of M), and
%   returns the new map. If WC is empty, you will be asked to select it
%   graphically. 
%   
%   See also DISKMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: center.m 196 2002-09-10 19:12:38Z driscoll $

if nargin == 1
  % Return center
  wc = map.center;
  if isempty(wc)
    p = polygon(map);
    wc = rsmap(0,vertex(p),angle(p)-1,...
        map.prevertex,map.prebranch,map.constant,map.qdata);
  end
  
else
  % Set center
  error('Recentering of RIESURF maps not supported.')
end
