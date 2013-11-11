function map = moebius(varargin)
%MOEBIUS Moebius transformation.
%   MOEBIUS(Z,W) creates the Moebius transformation the maps the
%   3-vector Z to W. Infinity is allowed in Z and W.
%   
%   MOEBIUS(C1,C2,C3,C4) creates the transformation
%   
%         C1 + C2*z
%         ---------
%         C3 + C4*z
%         
%   MOEBIUS([C1 C2 C3 C4]) also works.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: moebius.m 44 1998-07-01 20:21:54Z tad $

superiorto('double');

A = NaN*ones(1,4);

switch nargin
  case {0}
    map.source = [];
    map.image = [];
    map.coeff = [];
  case {1}
    C = varargin{1};
    if isa(C,'double') & length(C)==4
      map.source = [];
      map.image = [];
      map.coeff = C(:).';
    elseif isa(C,'moebius')
      map = C;
      return
    end
  case {4}
    map.source = [];
    map.image = [];
    map.coeff = cat(2,varargin{1:4});
  case {2}
    [z,w] = deal(varargin{1:2});
    % Pretty straightforward, but the infinities require care.
    if any(isinf(w))
      % Renumber to make w(2)=Inf
      j = find(isinf(w));
      renum = rem((j-1:j+1)+2,3)+1;
      z = z(renum);
      w = w(renum);
      % Depend on infinities in z
      if ~any(isinf(z))
        t1 = diff(z(1:2));
        t2 = -diff(z(2:3));
      elseif isinf(z(2))
        t1 = 1;
        t2 = 1;
      else
        if find(isinf(z))~=1
          % Move Inf to the beginning of z
          z = z([3 2 1]);
          w = w([3 2 1]);
        end
        A(1) = w(3)*(z(3)-z(2)) - w(1)*z(3);
        A(2) = w(1);
        A(3) = -z(2);
        A(4) = 1;
      end

    elseif any(isinf(z))
      % We already know all w are finite
      % Make z(2)=Inf
      j = find(isinf(z));
      renum = rem((j-1:j+1)+2,3)+1;
      z = z(renum);
      w = w(renum);
      t1 = -diff(w(2:3));
      t2 = diff(w(1:2));

    else					% everything finite
      t1 = -diff(z(1:2))*diff(w(2:3));
      t2 = -diff(z(2:3))*diff(w(1:2));
    end

    if isnan(A(1))
      % Did not encounter the special case above
      A(1) = w(3)*z(1)*t2 - w(1)*z(3)*t1;
      A(2) = w(1)*t1 - w(3)*t2;
      A(3) = z(1)*t2 - z(3)*t1;
      A(4) = t1 - t2;
    end

    map.source = z;
    map.image = w;
    map.coeff = A;
end

map = class(map,'moebius');

