function [w,beta,z,corners,renum] = rcorners(w,beta,z)
%RCORNERS (not intended for calling directly by the user)
%   Find corners of rectangle whose map is represented by prevertices z
%   on the strip, then renumber w, beta, and z (and the corners) so that
%   corners(1)=1.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rcorners.m 298 2009-09-15 14:36:37Z driscoll $

n = length(w);

% Deduce corner locations
left = abs(real(z)-min(real(z))) < eps;
right = abs(real(z)-max(real(z))) < eps;
top = abs(imag(z)-max(imag(z))) < eps;
bot = abs(imag(z)-min(imag(z))) < eps;
corners = find(left+right+top+bot - 1);
c1 = find(abs(z-max(real(z))) < eps);
offset = find(corners==c1);
corners = corners([offset:4,1:offset-1]);

% Renumber vertices so that corners(1)=1
renum = [corners(1):n,1:corners(1)-1];
w = w(renum);
beta = beta(renum);
z = z(renum);
corners = rem(corners-corners(1)+1+n-1,n)+1;

