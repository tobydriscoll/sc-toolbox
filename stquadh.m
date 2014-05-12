function I = stquadh(z1,z2,sing1,z,beta,qdat)

%   Copyright 1998 by Toby Driscoll.
%   $Id: stquadh.m 298 2009-09-15 14:36:37Z driscoll $

%   STQUAD applies the "1/2 rule" by assuming that the distance from the
%   integration interval to the nearest singularity is equal to the
%   distance from the left endpoint to the nearest singularity. This is
%   certainly true e.g. if one begins integration at the nearest
%   singularity to the target point. However, it may be violated in
%   important circumstances, such as when one integrates from the
%   next-nearest singularity (because the nearest maps to inf), or when
%   one is integrating between adjacent prevertices (as in the param
%   problem). The main difficulty is the possibility of singularities
%   from the "other side" of the strip intruding. 
%   
%   Here we assume that the integration intervals are horizontal. This
%   function recursively subdivides the interval until the "1/2 rule" is
%   satisfied for singularities "between" the endpoints. Actually, we
%   use a more stringent "alpha rule", for alpha > 1/2, as this seems to
%   be necessary sometimes.
% 
%   There must be no singularities *inside* the interval, of course.

n = length(z);
if isempty(sing1)
  sing1 = zeros(length(z1),1);
end

I = zeros(size(z1));
nontriv = find(z1(:)~=z2(:))';

for k = nontriv
  za = z1(k);
  zb = z2(k);
  sng = sing1(k);
 
  % alf==1/2 means the "1/2 rule." Better to be more strict.
  alf = .75;
  
  % Given integration length
  d = real(zb) - real(za);
  
  % Compute horizontal position (signed) and vertical distance (positive)
  % from the singularities to the left endpoint. If we are going from right
  % to left, reverse the sense of horizontal.
  dx = (real(z) - real(za)) * sign(d);
  dy = abs(imag(z) - imag(za));
  
  % We have to be concerned with singularities lying to the right (left if
  % d<0) of the left integration endpoint. 
  toright = (dx > 0) & (~isinf(z));
  % For points with small enough dx, the limitation is purely due to dy. For
  % others, it must be calculated.
  active = ( dx > dy/alf ) & toright;
  % Make sure that the left endpoint won't be included
  if sng
    active(sng) = 0;
  end
  
  % For those active, find the integration length constraint. This comes
  % from making the sing/right-endpoint distance equal to alf*L. 
  x = dx(active);
  y = dy(active);
  L = (x - sqrt((alf*x).^2 - (1-alf^2)*y.^2)) / (1-alf^2);
  
  % What is the maximum allowable integration length?
  L = min([ L(:); dy(toright & ~active)/alf ]);
  
  if L < abs(d)
    % Apply STQUAD on the safe part and recurse on the rest
    zmid = za + L*sign(d);
    I(k) = stquad(za,zmid,sng,z,beta,qdat);
    I(k) = I(k) + stquadh(zmid,zb,0,z,beta,qdat);
  else
    % No restriction
    I(k) = stquad(za,zb,sng,z,beta,qdat);
  end
  
end
  
