function [w,beta,aux] = scfix(type,w,beta,aux)
%SCFIX  Fix polygon to meet Schwarz-Christoffel toolbox constraints.
%   [W,BETA] = SCFIX(MAP,W,BETA) attempts to fix a problem in the given
%   polygon that arises from the posing of the parameter problem. SCFIX
%   is used when a call to xxPARAM results in an error and so
%   advises. In this case the polygon as given violates some fairly
%   arbitrary constraint. SCFIX remedies the situation by renumbering
%   the vertices, or, if necessary, adding a trivial (zero-turn)
%   vertex. MAP is one of {'hp','d','de','st','r'}. If one additional
%   input and output argument is given, it represents the indices of the
%   strip ends or the rectangle corners.
%   
%   See also SCCHECK, SCADDVTX.

%   Copyright 1998 by Toby Driscoll.
%   $Id: scfix.m 298 2009-09-15 14:36:37Z driscoll $

%   You may wonder, why not let the xxPARAM functions call SCFIX
%   automatically?  The trouble with that approach is that since a
%   function can't modify its inputs in the calling workspace,
%   either the xxPARAM functions would have to return more
%   arguments, or the mapping and plotting functions also would have
%   to detect and correct the problem every time they're called.
%   The problem is rare enough that this method seems adequate.

w = w(:);
beta = beta(:);
n = length(w);
renum = 1:n;

% Orientation conventions
sumb = -2 + 4*strcmp(type,'de');
if abs(sum(beta)+sumb) < 1e-9
  % Reverse order
  w = w([1,n:-1:2]);
  beta = scangle(w);
  renum = renum([1,n:-1:2]);
  if nargin > 3
    aux = renum(aux);
  end
end

% Less obvious restrictions
if strcmp(type,'hp') | strcmp(type,'d')
  shift = [2:n,1];
  % Renumber, if necessary, to meet requirements:
  %   w([1,2,n-1]) finite & sides at w(n) not collinear
  while any(isinf(w([1,2,n-1]))) | any(abs(beta(n)-[0,1])<eps)
    renum = renum(shift);
    w = w(shift);
    beta = beta(shift);
    if renum(1)==1 			% tried all orderings
      % First, be sure beta(n) is no longer a problem.
      if all((abs(beta-1)<eps)|(abs(beta)<eps))
	error('Polygon has empty interior!')
      end
      while any(abs(beta(n)-[0,1])<eps)
	w = w(shift);
	beta = beta(shift);
      end
      % Next, add one or two vertices as needed.
      if any(isinf(w(1:2)))
	[w,beta] = scaddvtx(w,beta,1);
	n = n+1;
      end
      if isinf(w(n-1))
	[w,beta] = scaddvtx(w,beta,n-1);
	n = n+1;
      end
      renum = 1:n;
      fprintf('\nWarning: A vertex has been added.\n\n')
      break
    end
  end
elseif strcmp(type,'de')
  shift = [2:n,1];
  % Renumber, if necessary, to ensure sides at w(n) not collinear
  %   (except if n==2, which is handled explicitly anyway)
  while any(abs(beta(n)-[0,1])<eps) & (n > 2)
    renum = renum(shift);
    w = w(shift);
    beta = beta(shift);
    if renum(1)==1
      deg = abs(beta) < eps;
      w(deg) = [];
      beta(deg) = [];
      renum = 1:2;
      n = 2;
      fprintf('\nPolygon is a line segment; removing superfluous vertices\n\n')
      break
    end
  end
elseif strcmp(type,'st')
  ends = aux;
  if isempty(ends)
    disp('Use mouse to select images of left and right ends of the strip.')
    figure(gcf)
    ends = scselect(w,beta,2);
  end
  renum = [ends(1):n,1:ends(1)-1];
  w = w(renum);
  beta = beta(renum);
  k = find(renum==ends(2));
  if k < 4
    if k < n-1
      % Switch ends.
      renum = [k:n,1:k-1];
      w = w(renum);
      beta = beta(renum);
      k = find(renum==1);
    else
      % Add one or two vertices.
      for j=1:4-k
	[w,beta] = scaddvtx(w,beta,j);
	n = n+1;
	k = k+1;
	fprintf('\nWarning: A vertex has been added.\n\n')
      end
    end
  end
  
  if k==n
    % Must add a vertex in any case.
    [w,beta] = scaddvtx(w,beta,n);
    n = n+1;
    fprintf('\nWarning: A vertex has been added.\n\n')
  end      

  
  if isinf(w(2))
    % Add two vertices.
    for j=1:2
      [w,beta] = scaddvtx(w,beta,j);
      n = n+1;
      k = k+1;
      fprintf('\nWarning: A vertex has been added.\n\n')
    end
  elseif isinf(w(3))
    % Add one vertex.
    [w,beta] = scaddvtx(w,beta,2);
    n = n+1;
    k = k+1;
    fprintf('\nWarning: A vertex has been added.\n\n')
  elseif isinf(w(n))
    [w,beta] = scaddvtx(w,beta,n);
    n = n+1;
    fprintf('\nWarning: A vertex has been added.\n\n')
  end
  
  aux = [1,k];

elseif strcmp(type,'r')
  corner = aux;
  renum = [corner(1):n,1:corner(1)-1];
  w = w(renum);
  beta = beta(renum);
  corner = rem(corner-corner(1)+1+n-1,n)+1;
  % Note: These problems are pretty rare.
  if any(abs(beta(n)-[0,1])<eps)
    % Try swapping sides 1-2 and 3-4.
    if ~any(abs(beta(corner(3)-1)-[0,1])<eps) & ~isinf(w(corner(3)))
      renum = [corner(3):n,1:corner(3)-1];
      w = w(renum);
      beta = beta(renum);
      corner = sort(rem(corner-corner(3)+1+n-1,n)+1);
    else
      error('Collinear sides make posing problem impossible')
    end
  end
  if isinf(w(2))
    [w,beta] = scaddvtx(w,beta,1);
    n = n+1;
    corner(2:4) = corner(2:4)+1;
  end

  aux = corner;
end
