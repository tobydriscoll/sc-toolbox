% Schwarz-Christoffel Toolbox
% Version 2.3   January 15, 2003
% Copyright (c) 1994-2003 by Toby Driscoll (driscoll@math.udel.edu).
%
% This is a selective listing of functions.
% See the user's guide for full usage details.
%
%   scdemo     - Select demos from a menu.
%   scgui      - Activate graphical user interface.
% 
% Polygons.
%   polygon    - Create a polygon object from vertices and angles.
%   polyedit   - Draw or edit a polygon with the mouse.
%   plot, fill - Plot a polygon, plot with filled interior.
%   display    - Show vertices.
%   length     - Number of vertices.
%   vertex     - Retrieve vertices.
%   angle      - Retrieve normalized angles.
%   (...)      - Reference or assign one or more vertices.
%   +, -, *    - Translate or scale.
%   winding    - Winding number around point(s).
%   linspace   - Equispace points around the boundary.
%   truncate   - Truncate infinite sides.
%   intersect  - Intersections of sides with a segment.
%   triangulate- Create a triangulation of the interior.
% 
% SC map types.
%   diskmap    - Disk -> polygon.
%   extermap   - Disk -> polygon exterior.
%   hplmap     - Half-plane -> polygon.
%   stripmap   - Strip -> polygon.
%   rectmap    - Rectangle -> generalized quadrilateral.
%   riesurfmap - Disk -> Riemann surface (multiple sheets).
%   crdiskmap  - Disk -> polygon, in cross-ratio formulation.
%   crrectmap  - Rectilinear polygon -> polygon, using cross-ratios.
% 
% Map operations.
%   eval, feval- Evaluate at point(s).
%   (...)      - Evaluate at point(s).
%   evalinv    - Evaluate inverse at point(s).
%   evaldiff   - Evaluate derivative at point(s).
%   plot       - Plot the image of an orthogonal grid under the map.
%   get        - Retrieve map parameters.
%   center     - Retrieve or set conformal center (disk maps only).
%   display    - Pretty-print.
%   polygon    - Retrieve the target polygon.
%   accuracy   - Approximate max-norm accuracy.
% 
% Additional utilities and applications.
%   scmapopt   - Set options for the parameter problem solution.
%   lapsolve   - Solve Laplace's equation using SC maps.
%   faber      - Compute Faber polynomials (of a polygon or extermap).
%   moebius    - Define a Moebius transformation.
%   composite  - Define a composition of maps.

