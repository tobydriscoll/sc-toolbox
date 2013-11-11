function [neww,orig,edge,triedge,edgetri] = crsplit(w)
%CRSPLIT Split polygon edges to ensure good crossratios.
%   CRSPLIT(W) splits the edges of polygon W to attempt to ensure that
%   subsequent triangulation will produce quadrilaterals having
%   reasonable crossratios. The first step is to "chop off" very sharp
%   exterior corners. Then, an iterative procedure subdivides the long
%   edges of narrow "channels." This process is guaranteed to terminate,
%   but not necessarily to be optimal.
% 
%   [WNEW,ORIG] = CRSPLIT(W) also returns a logical vector the same
%   length as WNEW. The 1 entries of ORIG correspond to the original
%   entries of W; i.e., WNEW(ORIG) returns the original W.
%       
%   CRSPLIT uses Delaunay triangulations to do its job.  [WNEW,E,TE,ET]
%   or [WNEW,ORIG,E,TE,ET] = CRSPLIT(...) returns a CDT as computed by
%   CRTRIANG/CRCDT.
%       
%   See also CRTRIANG, CRCDT.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crsplit.m 7 1998-05-10 04:37:19Z tad $

% First, isolate sharp corners
n = length(w);
beta = scangle(w);
sharp = (beta < -.75);
neww = NaN*zeros(3,n);
neww(2,:) = w(:).';
for j=find(sharp(:)')
  % Find distance to nearest "visible" vertex.
  jm1 = rem(j-2+n,n)+1;
  jp1 = rem(j,n)+1;
  v = (w-w(j))/(w(jp1)-w(j));
  ang = rem(angle(v)+2*pi,2*pi);
  vis = (ang>0) & (ang<ang(jm1));
  vis(rem([j-2:j]+n,n)+1) = [1;0;1];
  mindist = min(abs(w(vis)-w(j)));
  neww(1,j) = w(j) + 0.5*mindist*sign(w(jm1)-w(j));
  neww(3,j) = w(j) + 0.5*mindist*sign(w(jp1)-w(j));
end
w = neww(:);
sharp = [zeros(1,n);sharp(:)';zeros(1,n)];
sharp = sharp(:);
sharp(isnan(w)) = [];
orig = [zeros(1,n);ones(1,n);zeros(1,n)];
orig = orig(:);
orig(isnan(w)) = [];
w(isnan(w)) = [];
beta = scangle(w);

% Triangulate polygon. This is done to allow the use of geodesic distance.
[edge,triedge,edgetri] = crtriang(w);
[edge,triedge,edgetri] = crcdt(w,edge,triedge,edgetri);
% Build an adjacency matrix for the vertices in the triangulation
V = zeros(max(edge(:)));
for j=1:size(edge,2)
  V(edge(1,j),edge(2,j)) = 1;
  V(edge(2,j),edge(1,j)) = 1;
end

% Begin iterative phase
done = 0;
while ~done
  done = 1;
  n = length(w);
  neww = NaN*zeros(3,n);
  neww(1,:) = w(:).';
  for j=find(~sharp' & ~sharp([2:n,1])')
    % Vertices in sequence, mod n
    jp1 = rem(j,n)+1;
    jm1 = rem(j-2+n,n)+1;
    jp2 = rem(jp1,n)+1;
    
    % All points other than j and jp1
    mask = ones(n,1);
    mask([j,jp1]) = zeros(2,1);
    % Points adjacent to j in triangulation
    mask1 = V(:,j) & mask;
    % Exclusion in the case of a slit (would appear as zero distance)
    if abs(beta(j)-1) < eps
      mask1(jm1) = 0;
    end
    w1 = w(mask1);
    % Points adjacent to jp1, with slit exclusion
    mask2 = V(:,jp1) & mask;
    if abs(beta(jp1)-1) < eps
      mask2(jp2) = 0;
    end
    w2 = w(mask2);
    % Geodesic distance from adjacent points to segment [w(j),w(jp1)]
    [dist1,dist2] = crpsgd(w([j,jp1]),w1,w2);
    dist = [dist1(:);dist2(:)];

    % If distance is too small, subdivide
    if min(dist)/abs(w(jp1)-w(j)) < 1/(sqrt(2)*3)
      done = 0;
      neww(2:3,j) = w(j) + [1/3;2/3]*(w(jp1)-w(j));
    end
  end
  
  % Update triangulation
  % First, update edge vertex numbers
  renum = neww;
  renum(~isnan(neww)) = 1:sum(~isnan(neww(:)));
  edge(:) = renum(1,edge);
  % Now update triangles & edges
  new = find(~isnan(neww(2,:)));
  newn = sum(~isnan(neww(:)));
  for j=new
    idx = renum(1,[j rem(j,n)+1]);	% endpoints of newly split edge
    % Find split edge
    e = find((edge(1,:)==idx(1)) & (edge(2,:)==idx(2)));
    if isempty(e),
      e = find((edge(1,:)==idx(2)) & (edge(2,:)==idx(1)));
    end 
    t = edgetri(1,e);			% triangle that edge e was in
    edge(:,e) = idx(1)+[0;1];		% Replace e by new edge
    te = triedge(:,t);			% edges of t
    m = find(te==e);			% e's place among t's edges
    ntri = size(triedge,2);
    nedge = size(edge,2);
    % Find the edge in triangle t that is ccw from e
    i2 = rem(m,3)+1;
    e2 = te(i2);
    if ~any(edge(:,e2)==idx(2))
      i2 = rem(m+1,3)+1;
      e2 = te(i2);
    end
    vtx = edge(edge(:,e2)~=idx(2),e2);	% 3rd vertex of t
    % New interior edges
    edge(:,nedge+(1:2)) = [[idx(1)+1;vtx] [idx(1)+2;vtx]];
    % New boundary edges
    edge(:,nedge+(3:4)) = rem(idx(1)+[[1;2] [2;3]]-1,newn)+1;
    % 2 new triangles and an old one replaced
    triedge(:,ntri+(1:2)) = [[nedge+[3;2;1]] [nedge+2;nedge+4;e2]];
    triedge(i2,t) = nedge+1;
    % New triangle memberships for new edges and e2
    edgetri(:,nedge+(1:2)) = [[t;ntri+1] [ntri+1;ntri+2]];
    edgetri(:,nedge+(3:4)) = [[ntri+1;0] [ntri+2;0]];
    edgetri(edgetri(:,e2)==t,e2) = ntri+2;
  end
  
  w = neww(:);
  sharp = [sharp(:)';zeros(2,n)];
  sharp = sharp(:);
  sharp(isnan(w)) = [];
  orig = [orig(:)';zeros(2,n)];
  orig = orig(:);
  orig(isnan(w)) = [];
  w(isnan(w)) = [];
  beta = scangle(w);

  % Recompute CDT for geodesic distance
  [edge,triedge,edgetri] = crcdt(w,edge,triedge,edgetri);
  % Build an adjacency matrix for the vertices in the triangulation
  V = zeros(max(edge(:)));
  for j=1:size(edge,2)
    V(edge(1,j),edge(2,j)) = 1;
    V(edge(2,j),edge(1,j)) = 1;
  end
  
end

orig = logical(orig);

neww = w;
% Special output format
if nargout == 4
  orig = edge;
  edge = triedge;
  triedge = edgetri;
end


function [d1,d2] = crpsgd(segment,pts1,pts2)
%CRPSGD Geodesic distance from points to a line segment.

%   [D1,D2] = CRPSGD(SEG,PTS1,PTS2) computes geodesic distance from two
%   sets of points, PTS1 and PTS2, to a directed segment SEG. "Geodesic
%   distance" here means that the segment describes an "inside" (to the
%   left along the directed segment) and an "outside". For points
%   "inside", the distance is the usual point-segment distance. For
%   points "outside", the distance for a point in PTS{j} is the distance
%   to SEG(j). Typically the points in PTS{j} are vertices of a polygon,
%   adjacent to SEG(j) in a constrained Delaunay triangulation, and one
%   needs to know the polygon-interior distance from each vetex to the
%   polygon side SEG.

d1 = Inf*pts1;
d2 = Inf*pts2;

rot = 1/sign(diff(segment));

% For points to the left of the segment, use normal pt-seg distance. For
% others, use distance to the vertex it "belongs" to.

pts = (pts1-segment(1))*rot;
mask = imag(pts) > eps;
if sum(mask)
  d1(mask) = crpsdist(segment,pts1(mask));
end
d1(~mask) = abs(pts1(~mask)-segment(1));

pts = (pts2-segment(1))*rot;
mask = imag(pts) > eps;
if sum(mask)
  d2(mask) = crpsdist(segment,pts2(mask));
end
d2(~mask) = abs(pts2(~mask)-segment(2));

