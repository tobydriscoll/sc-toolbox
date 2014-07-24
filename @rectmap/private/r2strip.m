function [yp,yprime] = r2strip(zp,z,L)
%R2STRIP Map from rectangle to strip.
%   R2STRIP(ZP,Z,L) maps from a rectangle to the strip 0 <= Im z <= 1,
%   with the function log(sn(z|m))/pi, where sn is a Jacobi elliptic
%   function and m = exp(-2*pi*L).  The prevertices of the map (in the
%   rectangle domain) are given by Z; only the corners of the rectangle
%   defined by Z are used.
% 
%   The derivative of the map is returned as a second argument.
%
%   NOTE: The functionality is NOT parallel to HP2DISK and DISK2HP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: r2strip.m 298 2009-09-15 14:36:37Z driscoll $

% The built-in ellipj accepts only real arguments. The standard identity is
% not helpful when m is near zero, because (1-m) loses digits of accuracy.

K = max(real(z));
Kp = max(imag(z));
yp = zp;
yprime = zp;

[sn,cn,dn] = ellipjc(zp,L);
% Make sure everything is in the upper half-plane (fix roundoff)
sn = real(sn) + i*max(imag(sn),0);

yp = log(sn)/pi;
yprime = cn.*dn./sn/pi;

% Make sure everything is in the strip (roundoff could put it outside)
yp = real(yp) + i*max(0,imag(yp));
yp = real(yp) + i*min(1,imag(yp));

