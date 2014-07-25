function [knz,inz] = nearz(z,dataz)
% NEARZ determines the nearest vertex to a given pt z. On return integer
% inz indicates that the nearest pt is found either in Z0 (inz = 0) or
% in Z1 (inz = 1). knz contains the corresponding index in Z0 or in Z1.
% If on return inz=2, that means no appropriate vertex is found, which is a
% device used for inverse mapping routine.
inz = 2;
dist = 99;

d = abs(z - dataz.Z0);
knz = find(dataz.ALFA0 > 0 & d < dist);
if isempty(knz)==0
    knz = find(d==min(d(knz)));
    knz = knz(1); % Avoid multiple nearest vertices.
    dist = d(knz);
    inz = 0;
end

d = abs(z - dataz.Z1);
idx = find(d < dist);
if isempty(idx)==0
    knz = find(d==min(d(idx)));
    knz = knz(1); % Avoid multiple nearest vertices.
    inz = 1;
end
