function [A,B,C,D] = moebius(z,w)
%MOEBIUS Moebius transformation parameters.
%   A = MOEBIUS(Z,W) computes the coefficients of the Moebius
%   transformation taking the 3-vector Z to W, so that
%   
%            W = (A(1)*Z + A(2))./(A(3)*Z + A(4)).
%   
%   Infinities are allowed.
% 
%   If four output arguments are used, they will be given the A values.

%   Copyright 1998 by Toby Driscoll.
%   $Id: moebius.m 298 2009-09-15 14:36:37Z driscoll $

A = NaN*ones(1,4);

if any(isinf(w))
  % Make w(2)=Inf
  j = find(isinf(w));
  renum = rem((j-1:j+1)+2,3)+1;
  z = z(renum);
  w = w(renum);
  if ~any(isinf(z))
    t1 = diff(z(1:2));
    t2 = -diff(z(2:3));
  elseif isinf(z(2))
    t1 = 1;
    t2 = 1;
  else
    % Will have to deal separately with w(2)=Inf and z(1)=Inf
    j = find(isinf(z));
    if j~=1
      z = z([3 2 1]);
      w = w([3 2 1]);
    end
    A(1) = w(1);
    A(2) = w(3)*(z(3)-z(2)) - w(1)*z(3);
    A(3) = 1;
    A(4) = -z(2);
  end

elseif any(isinf(z))
  % We already know ~any(isinf(w))
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
  A(1) = w(1)*t1 - w(3)*t2;
  A(2) = w(3)*z(1)*t2 - w(1)*z(3)*t1;
  A(3) = t1 - t2;
  A(4) = z(1)*t2 - z(3)*t1;
end

if nargout==4
  D = A(4);
  C = A(3);
  B = A(2);
  A = A(1);
end
