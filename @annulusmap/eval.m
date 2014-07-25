function z = eval(map,w)
%EVAL Evaluate Schwarz-Christoffel annulus map at points.
%   EVAL(MAP,Z) evaluates the Schwarz-Christoffel map MAP at the points Z
%   in the canoncial annulus. 
%
%   See also ANNULUSMAP.

%   Copyright 2014 by Toby Driscoll.
%   Written by Alfa Heryudono.

% TODO: check if w is in the annulus

kww = 0;
ic = 2;

% check if w is in W0.
idx = find(map.w0 == w, 1 );
if isempty(idx)==0
    kww = idx;
    ic = 0;
else
    % check if w is in W1.
    idx = find(map.w1 == w, 1 );
    if isempty(idx)==0
        kww = idx;
        ic = 1;
    end
end
nptq = 8;

%Making the bridge to old subroutine
dataz = struct('M',map.M,'N',map.N,'Z0',map.Z0,'Z1',map.Z1,'ALFA0',map.ALFA0,'ALFA1',map.ALFA1,'ISHAPE',map.ISHAPE);
z = zdsc(w,kww,ic,map.u,map.c,map.w0,map.w1,map.phi0,map.phi1,nptq,map.qwork,1,dataz);

end