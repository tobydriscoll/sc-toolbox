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
%   $Id: center.m 30 1998-06-22 22:46:23Z tad $

if nargin == 1
  % Return center
  wc = map.center;
  if isempty(wc)
    p = polygon(map);
    wc = dmap(0,vertex(p),angle(p)-1,...
        map.prevertex,map.constant,map.qdata);
  end
  
else
  % Set center
  p = polygon(map);
  qdata = map.qdata;
  z = map.prevertex;
  w = vertex(p);
  beta = angle(p) - 1;
  
  if isempty(wc)
    fig = figure;
    plot(p)
    title('Click at conformal center')
    [xc,yc] = ginput(1);
    wc = xc + i*yc;
    delete(fig)
  end

  if ~any(isinf(w)) & ~isinpoly(wc,p)
    error('Conformal center must be inside polygon.')
  end
  
  % Find inverse image of wc under current map
  zc = dinvmap(wc,w,beta,z,map.constant,qdata);

  % Use Moebius transform to reset prevertices
  y = ((1-zc')/(1-zc))*(z-zc)./(1-zc'*z);
  y(length(y)) = 1;			% force it to be exact
  y = y./abs(y);
  
  % Recalculate constant
  mid = mean(y(1:2));
  I = dquad(y(1),mid,1,y,beta,qdata) - dquad(y(2),mid,2,y,beta,qdata);
  c = diff(w(1:2))/I;
  
  % Assign new values
  map.prevertex = y;
  map.constant = c;
  map.center = wc;
  map.accuracy = [];
  map.accuracy = accuracy(map);
  wc = map;
end
