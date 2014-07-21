function err = sccheck(type,w,beta,aux)
%SCCHECK Check polygon inputs to Schwarz-Christoffel functions.
%   SCCHECK(MAP,W,BETA) is used by the xxPARAM functions to check the
%   validity of inputs describing the polygon to be mapped.  MAP is a
%   string consisting of the prefix to PARAM ('d', 'dp', etc.).
%   
%   If errors are found, execution will terminate.  Sometimes the
%   trouble has to do with how the parameter problem is posed, which
%   imposes a few nonobvious constraints.  The function SCFIX is
%   provided to automatically fix such difficulties, by renumbering or
%   perhaps adding vertices. SCCHECK output is 1 if the problem is
%   rectifiable by SCFIX, 2 if warning status only.
%   
%   See also SCFIX.

%   Copyright 1998 by Toby Driscoll.
%   $Id: sccheck.m 298 2009-09-15 14:36:37Z driscoll $

w = w(:);
beta = beta(:);
n = length(w);
atinf = logical(beta <= -1);
renum = 1:n;
err = 0;

% Universal truths
if length(beta)~=n
  error('Mismatched angles and vertices')
elseif any(beta > 1) | any(beta < -3)
  error('Each angle must be in [-3,1]')
end

% Infinite vertices
if ~strcmp(type,'de') & ~strcmp(type,'cr')
  if any(isinf(w(~atinf))) | any(~isinf(w(atinf)))
    error('Infinite vertices must correspond to angle <= -1')
  else
    da = diff(find(atinf));
    if ~isempty(da) & any(da==1)
      error('Infinite vertices must not be adjacent')
    end
  end
else
  if any(atinf) | any(isinf(w))
    error('Infinite vertices not allowed')
  end
end
sumb = -2 * (-1)^(strcmp(type,'de'));

% Orientation conventions
if abs(sum(beta)+sumb) < 1e-9
  fprintf('\nVertices were probably given in wrong order\n')
  err = 1;
elseif abs(sum(beta)-sumb) > 1e-9
  fprintf('\nWarning: Angles do not sum to %d\n\n',sumb)
  err = 2;
end

% Some finer points
if strcmp(type,'hp') | strcmp(type,'d')
  if n < 3
    error('Polygon must have at least three vertices')
  elseif any(isinf(w([1,2,n-1])))
    fprintf('\nInfinite vertices must not be at positions 1, 2, or n-1\n')
    err = 1;
  elseif any(abs(beta(n)-[0,1])<eps)
    fprintf('\nSides adjacent to w(n) must not be collinear\n')
    err = 1;
  end
elseif strcmp(type,'cr')
  if n < 4
    error('Polygon must have at least four vertices')
  end
elseif strcmp(type,'de')
  if n < 2
    error('Polygon must have at least two vertices')
  elseif (beta(n)==0 | beta(n)==1) & (n > 2)
    fprintf('\nSides adjacent to w(n) must not be collinear\n')
    err = 1;
  end
elseif strcmp(type,'st')
  if n < 5
    error('Polygon must have at least five vertices')
  end
  ends = aux;
  renum = [ends(1):n,1:ends(1)-1];
  w = w(renum);
  beta = beta(renum);
  k = find(renum==ends(2));
  if any(atinf([2,3,n]))
    fprintf('\nVertices at (w(ends(1)) + [1,2,-1]) must be finite\n')
    err = 1;
  elseif k-2 < 2
    fprintf('\nThere must be at least 2 vertices between ends 1 and 2\n')
    err = 1;
  elseif k==n
    fprintf('\nThere must be at least one vertex between ends 2 and 1\n')
    err = 1;
  end
elseif strcmp(type,'r')
  corner = aux;
  renum = [corner(1):n,1:corner(1)-1];
  w = w(renum);
  beta = beta(renum);
  corner = rem(corner-corner(1)+1+n-1,n)+1;
  if n < 4
    error('Polygon must have at least four vertices')
  elseif corner~=sort(corner)
    error('Corners must be specified in ccw order')
  elseif isinf(w(1))
    error('Corner(1) must be finite')
  end
  if isinf(w(2))
    fprintf('\nVertex corner(1)+1 must be finite\n')
    err = 1;
  end
  if any(abs(beta(n)-[0,1])<eps)
    fprintf('\nSides adjacent to w(corner(1)-1) must not be collinear\n')
    err = 1;
  end
   
end
