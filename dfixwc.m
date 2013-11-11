function [y,d] = dfixwc(w,beta,z,c,wc,tol)
%DFIXWC Fix conformal center of disk map.
%   The conformal center WC of a Schwarz-Christoffel interior disk map
%   is defined as the image of zero.  The parameter problem solver
%   DPARAM does not allow control over the placement of the conformal
%   center.  Using the output Z,C from DPARAM, [Z0,C0] =
%   DFIXWC(W,BETA,Z,C,WC) computes a Moebius transformation so that if
%   Z0 and C0 are used in place of Z and C, the conformal center of the
%   resulting map will be WC.
%   
%   [Z0,C0] = DFIXWC(W,BETA,Z,C,WC,TOL) uses tolerance TOL.
%
%   See also DPARAM, PTSOURCE.
  
%   Copyright 1998 by Toby Driscoll.
%   $Id: dfixwc.m 7 1998-05-10 04:37:19Z tad $

n = length(w);

if nargin < 6
  [trace,tol,method] = scparopt([]);
end 

zc = dinvmap(wc,w,beta,z,c,tol);

% Transform prevertices.
y = ((1-zc')/(1-zc))*(z-zc)./(1-zc'*z);
y(n) = 1;				% force it to be exact
y = y./abs(y);

% Recalculate constant from scratch.
mid = (y(1)+y(2))/2;
qdat = scqdata(beta,ceil(-log10(tol)));
d = (w(1) - w(2))/...
    (dquad(y(2),mid,2,y,beta,qdat) - dquad(y(1),mid,1,y,beta,qdat));

