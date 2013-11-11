function [wp,qnum,up] = crrmap(zp,w,beta,wr,betar,cr,aff,affr,Q,qdat,qdatr)
%CRRMAP Schwarz-Christoffel rectified map in crossratio formulation.
%   CRRMAP(ZP,W,BETA,WR,BETAR,CR,AFF,AFFR,Q,QDAT,QDATR) computes the
%   values of the rectified map from WR to W at the points in vector
%   ZP. The arguments are returned from CRPARAM, CRRECT, and SCQDATA.
%
%   If a scalar TOL is supplied in place of the last two arguments,
%   quadrature data intended to give an answer accurate to within
%   roughly TOL is used. Omitting these arguments altogether leads to
%   TOL = 1e-8.
% 
%   It is possible that the rectified polygon lies on multiple Riemann
%   sheets; i.e., overlaps itself. Hence the map may be multiply
%   valued. If at least one point of ZP lies in more than one covering,
%   then each row of output will contain all possible values for the map
%   at the corresponding point of ZP. Unused entries will be assigned
%   NaN. 
%   
%   [WP,QNUM,UP] = CRRMAP(...) also returns the indices of the
%   quadrilaterals used to compute the maps, and the intermediate points
%   found in the disk, in those embeddings.
%       
%   Note that by switching the roles of W, BETA, and AFF with WR, BETAR,
%   and AFFR, one inverts the map.
%
%   See also CRPARAM, CRRECT, CRRPLOT, CRRINVMP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crrmap.m 7 1998-05-10 04:37:19Z tad $

% Parse input and initialize
n = length(w);
w = w(:);
beta = beta(:);

if nargin < 11
  tol = 1e-8;
  qdat = scqdata(beta,8);
  qdatr = scqdata(betar,8);
elseif length(qdat)==1
  tol = qdat;
  qdat = scqdata(beta,max(ceil(-log10(tol)),3));
  qdatr = scqdata(betar,max(ceil(-log10(tol)),3));
else
  tol = 10^(-size(qdat,1));
end

% Look for multiple-sheet case
[zpindex,onvtx] = isinpoly(zp(:),wr,betar,tol);
p = length(zp(:));
maxindex = max(zpindex);
if maxindex == 0
  % zp consists only of points outside the domain
  wp = NaN*zp;
  return
elseif maxindex == 1
  shape = size(zp);
else
  shape = [p maxindex];
end
numcopies = zeros(p,1);		% # of copies found so far per point
zp = zp(:);
wp = NaN*zeros(p,maxindex);
up = wp;
qnum = NaN*zeros(p,maxindex);	% where each point was mapped from

% First, deal with points coincident with vertices
for k = find(any(onvtx))
  vtxnum = find(onvtx(:,k));
  nv = length(vtxnum);
  wp(k,1:nv) = w(vtxnum).';
  for j = 1:nv
    qn(k,j) = min(find( any(vtxnum(j)==Q.qlvert) ));
  end
  numcopies(k) = nv;
end  

% For each embedding, compose inverse rectified & forward original maps for
% points of zp in the quadrilateral
for q=1:n-3
  % Short-circuit if all is done
  if all(numcopies == zpindex), break, end

  % Still need to be mapped
  idx = find(numcopies < zpindex);

  % Those that are inside quad'l q
  mask = logical(abs(isinpoly(zp(idx),wr(Q.qlvert(:,q)),tol)));

  % If a point was mapped from a neighboring quad'l, exclude it, since it
  % would be a repeat
  if maxindex > 1
    qn = qnum(idx,:)';
    check = ~isnan(qn);
    qn(check) = Q.adjacent(q,qn(check));
    mask = mask & ~any(qn==1)';
  end
  
  if any(mask)
    % Proceed
    idx = idx(mask);
    z = crembed(cr,Q,q);
    up1 = crimap0(zp(idx),z,betar,affr(q,:),qdatr,[0 tol]);
    wp1 = crmap0(up1,z,beta,aff(q,:),qdat);
    % As a further check on repeats, compare values to previous ones
    % (Points may lie on diagonals and appear to be in non-neighboring
    % quadrilaterals)
    if maxindex > 1
      dif = abs(wp1(:,ones(maxindex,1)) - wp(idx,:));
      unique = ~any(dif < 10*tol,2);
      wp1 = wp1(unique);
      up1 = up1(unique);
      idx = idx(unique);
    end
    % Fill in new values
    new = idx + numcopies(idx)*p;
    wp(new) = wp1;
    up(new) = up1;
    qnum(new) = q*ones(length(idx),1);
    numcopies(idx) = numcopies(idx) + 1;
      
  end
end

wp = reshape(wp,shape(1),shape(2));
up = reshape(up,shape(1),shape(2));
