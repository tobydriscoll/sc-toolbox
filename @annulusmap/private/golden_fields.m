function s = golden_fields(map)
%GOLDEN_FIELDS (not intended for calling directly by the user)
%   Extracts annulusmap internal fields for the C++ golden-value test
%   generators (tests/cpp/generators), since annulusmap's custom subsref
%   blocks ordinary dot-notation property access from outside the class.
s = struct('M', map.M, 'N', map.N, 'ALFA0', map.ALFA0, 'ALFA1', map.ALFA1, ...
    'u', map.u, 'c', map.c, 'w0', map.w0, 'w1', map.w1, ...
    'phi0', map.phi0, 'phi1', map.phi1);
end
