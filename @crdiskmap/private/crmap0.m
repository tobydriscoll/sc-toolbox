function wp = crmap0(zp,z,beta,aff,qdat)
%CRMAP Single-embedding map in crossratio formulation.
%   CRMAP0(ZP,Z,BETA,AFF) computes the image of ZP under the map defined
%   by the single prevertex embedding Z and the affine transformation
%   AFF(1:2).
%       
%   CRMAP0(ZP,Z,BETA,AFF,TOL) uses quadrature data intended to give an
%   answer accurate to within roughly TOL.
%       
%   CRMAP0(ZP,Z,BETA,AFF,WC) uses a tolerance of 1e-8.
% 
%   In keeping with the CR approach, the integration is from the center
%   from the disk. Results may not be accurate for points near crowded
%   prevertices. Instead one should re-embed.
%
%   See also CRPARAM, CREMBED, CRAFFINE, CRMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crmap0.m 7 1998-05-10 04:37:19Z tad $

% Parse input and initialize
z = z(:);
n = length(z);
beta = beta(:);
if nargin < 5
  qdat = sctool.scqdata(beta,8);
elseif length(qdat)==1
  qdat = sctool.scqdata(beta,max(ceil(-log10(qdat)),2));
end
wp = zp;
zp = zp(:);
np = length(zp);

% Single out points that are essentially coincident with a prevertex
dif = abs(z(:,ones(np,1)) - zp(:,ones(n,1)).') < 10*eps;
[ir,ic] = find(dif);
sing = zeros(np,1);
% Assign them accurate Gauss-Jacobi quadrature
sing(ic) = ir(:);

% Do the maps
wp(:) = -aff(1)*crquad(zp,sing,z,beta,qdat) + aff(2);
