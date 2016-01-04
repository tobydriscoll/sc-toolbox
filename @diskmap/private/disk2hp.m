function [zhp,chp] = disk2hp(w,beta,z,c)
%DISK2HP Convert solution from the disk to one from the half-plane.
%   [ZHP,CHP] = DISK2HP(W,BETA,Z,C) quickly transforms the solution Z,C
%   of the Schwarz-Christoffel disk mapping parameter problem to the
%   solution ZHP,CHP of the half-plane problem.
%   
%   See also HP2DISK, DPARAM, HPPARAM.
      
%   Copyright 1998 by Toby Driscoll.
%   $Id: disk2hp.m 7 1998-05-10 04:37:19Z tad $

n = length(w);
zhp = zeros(size(z));
zhp(n) = Inf;
zhp(1:n-1) = -i*(z(1:n-1)+1)./(z(1:n-1)-1); % Mobius transfmn
zhp = real(zhp);	

% Recalculate constant from scratch.
mid = mean(zhp(1:2));
qdat = sctool.scqdata(beta(1:n-1),16);
chp = (w(1)-w(2))/(hplmap.hpquad(zhp(2),mid,2,zhp(1:n-1),beta(1:n-1),qdat) - ...
    hplmap.hpquad(zhp(1),mid,1,zhp(1:n-1),beta(1:n-1),qdat));


