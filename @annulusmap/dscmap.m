function z = dscmap(w,map)
% DSCMAP maps a pt w to a pt z

% DSCMAP(w,map) maps a pt w in the annulus to the doubly connected region 
% map is an annulusmap object.

% see @annulusmap/annulusmap.m

% check if w is in the annulus
% Not implemented yet.
%----------------------
%----------------------

kww = 0;
ic = 2;
% check if w is in W0.
idx = min(find(map.w0 == w));
if isempty(idx)==0
    kww = idx;
    ic = 0;
else
    % check if w is in W1.
    idx = min(find(map.w1 == w));
    if isempty(idx)==0
        kww = idx;
        ic = 1;
    end
end
nptq = 8;

%Making the bridge to old subroutine
dataz = struct('M',map.M,'N',map.N,'Z0',map.Z0,'Z1',map.Z1,'ALFA0',map.ALFA0,'ALFA1',map.ALFA1,'ISHAPE',map.ISHAPE);
z = zdsc(w,kww,ic,map.u,map.c,map.w0,map.w1,map.phi0,map.phi1,nptq,map.qwork,1,dataz);