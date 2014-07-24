function [zp,wp,new] = scpadapt(zp,wp,minlen,maxlen,clip)
%SCPADAPT (not intended for calling directly by the user)

%   This function looks for unacceptable nonsmoothness in the curve(s)
%   represented by wp. At such points it adds in-between points to
%   zp. On return zp is the combined vector of points in order, wp has
%   NaN's at the new points, and new is a 0-1 vector flagging the newly
%   added points.
%       
%   The algorithm is basically stolen from FPLOT. If extrapolation of
%   the linear interpolation that the graphics will use results in an
%   estimate too far from reality, refinement is called for.
% 
%   When the clip argument is used, the criteria do not apply for points
%   outside the clip region. Why refine for invisible points?
%       
%   wp can have multiple columns (curves). Any curve needing refinement
%   causes it for all of them.
%       
%   This should become the adaptive routine for all SC
%   plotting. allowances must be made for cases when the ends of zp need
%   not be bounded (as with the half-plane and strip). Perhaps the newly
%   added points should be spaced exponentially rather than
%   algebraically.
       
%   Copyright 1997-2003 by Toby Driscoll. 
%   $Id: scpadapt.m 252 2003-03-07 16:34:55Z driscoll $

m = size(wp,1);
dwp = diff(wp);

% Infinities at the ends of zp mean that we could go out
% forever. However, we will pick a finite bound.
linf = isinf(zp(1));
if linf
  zp(1) = zp(2) + max(3,abs(zp(3)-zp(2))^2)*sign(zp(2)-zp(3));
end
rinf = isinf(zp(m));
if rinf
  zp(m) = zp(m-1) + max(3,abs(zp(m-1)-zp(m-2))^2)*sign(zp(m-1)-zp(m-2));
end

% Sines of the turning angles at interior points
sines = abs(sin(angle(dwp(1:m-2,:).*conj(dwp(2:m-1,:))) ));

% Distances from linear extrapolation to actual value. Each interior point
% has a column, with two rows being errors to either side.
dwp = abs(dwp);
err = [ max(dwp(1:m-2,:).*sines,[],2) max(dwp(2:m-1,:).*sines,[],2) ];

% Forget about NaNs.
err(isnan(err)) = 0;
dwp(isnan(dwp)) = 0;

% Flag large errors
bad = logical([ 0 0; (err > maxlen/12); 0 0 ]);

% Also flag if the segments themselves are too long
% However, exclude where segments are very short
long = any(dwp > maxlen, 2);
short = all(dwp < minlen, 2);
bad(:,1) = (bad(:,1) | [0; long]) & [0;~short];
bad(:,2) = (bad(:,2) | [long; 0]) & [~short;0]; 


% Find clipped points 
xp = real(wp);
yp = imag(wp);
lohi = (yp < clip(3)) | (yp > clip(4)) | isinf(wp);
lfrt = (xp < clip(1)) | (xp > clip(2)) | isinf(wp);

% To be refined, you or a neighbor must be inside
inside = all(~lohi & ~lfrt,2);
bad(:,1) = bad(:,1) & (inside | [1;inside(1:m-1)]);
bad(:,2) = bad(:,2) & (inside | [inside(2:m);  1]);

% Exception: Extend the lines that bound the axes box. If you cross one
% horizontal and one vertical among those lines between adjacent outside
% points, you must refine in-between. Otherwise you could miss valid inside
% points, or get a line that crosses back into the axes box.
lh = lohi & ~lfrt;
lr = lfrt & ~lohi;
flag = (lh(1:m-1,:) & lr(2:m,:)) | (lr(1:m-1,:) & lh(2:m,:));
flag = any(flag,2);
if any(flag)
  bad(find(flag)+1,1) = 1;
  bad(find(flag),2) = 1;
end


% Refinement
nans = repmat(NaN,1,m);
zp = [nans; zp(:).'; nans];
new = [nans; zeros(1,m); nans];

% Add new source points. Always go one-third of the distance toward the
% offending party.
back = find(bad(:,1));
if ~isempty(back)
  zp(1,back) = (zp(2,back-1) + 2*zp(2,back)) / 3;
  new(1,back) = 1;
end
forw = find(bad(:,2));
if ~isempty(forw)
  zp(3,forw) = (zp(2,forw+1) + 2*zp(2,forw)) / 3;
  new(3,forw) = 1;
end
  
% If we had infinites at the ends, put them back
if linf
  zp(2,1) = Inf;
end
if rinf
  zp(2,m) = Inf;
end

% Restack zp into a vector, deleting NaNs
mask = ~isnan(zp);
zp = zp(mask);
new = new(mask);
  
% Update wp
wpold = wp;
wp = repmat(NaN,length(zp),size(wpold,2));
wp(~new,:) = wpold;

new = logical(new);
