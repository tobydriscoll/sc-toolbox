function [wcfix,wc] = crfixwc(w,beta,cr,aff,Q,wc)
%CRFIXWC Fix conformal center in crossratio formulation.
%   The conformal center is defined as the image of zero from the
%   disk. The parameter problem solution obtained from CRPARAM
%   deliberately leaves this unspecified. Before one can compute forward
%   or inverse maps, then, one must first fix the conformal center.
%       
%   WCFIX = CRFIXWC(W,BETA,CR,AFF,Q,WC) returns a vector containing the
%   information that CRMAP and CRINVMAP need in order to place the
%   conformal center at WC. You must first run CRAFFINE to find AFF.
%   
%   If you leave out the argument WC, you will be prompted for it
%   graphically, and it will be the second output argument.  
%       
%   See also CRPARAM, CRAFFINE, CRMAP, CRINVMAP.
 
%   Copyright 1998 by Toby Driscoll.
%   $Id: crfixwc.m 7 1998-05-10 04:37:19Z tad $

n = length(cr)+3;

if nargin < 6
  fig = figure;
  plotpoly(w,beta)
  title('Click to set conformal center')
  [x,y] = ginput(1);
  wc = x+i*y;
  delete(fig)
  drawnow
end

% Find a quadrilateral containing wc
for quadnum=1:n-3
  if sctool.isinpoly(wc,w(Q.qlvert(:,quadnum)),1e-4)
    break
  end
end

% Invert for that embedding
z = crembed(cr,Q,quadnum);
zc = crimap0(wc,z,beta,aff(quadnum,:));

% Find the Moebius transform that goes from the original disk to this
% embedding. Always assume that original 1 maps to w(n).
a = sign((1-zc'*z(n))/(z(n)-zc));
mt = [1 zc*a zc' a];

wcfix = [quadnum mt];
