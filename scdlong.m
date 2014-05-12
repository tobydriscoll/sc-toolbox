function slide=scdlong
% This is a slideshow file for use with playshow.m and makeshow.m
% To see it run, type 'playshow scdlong', 

% Copyright (c) 1984-98 by The MathWorks, Inc.
% $Id: scdlong.m 298 2009-09-15 14:36:37Z driscoll $
if nargout<1,
  playshow scdlong
else
  %========== Slide 1 ==========

  slide(1).code={
   '' };
  slide(1).text={
   'Press "Start" to see a demonstration of mapping to elongated regions.'};

  %========== Slide 2 ==========

  slide(2).code={
   'p = polygon([ 2.0000 + 3.0000i',
   '  -3.0000 + 3.0000i',
   '  -3.0000 - 3.0000i',
   '  -2.0000 - 2.0000i',
   '  -2.0000 + 2.0000i',
   '   2.0000 + 2.0000i  ]);',
   'plot(p)',
   'f=diskmap(p,scmapopt(''tol'',1e-6,''trace'',''on''));',
   '',
   '',
   '' };
  slide(2).text={
   'Here is a mildly elongated region.',
   '>> plot(p)',
   '',
   'We find its SC disk map.',
   '>> f = diskmap(p);'};

  %========== Slide 3 ==========

  slide(3).code={
   'z=prevertex(f);',
   'plot(exp(i*linspace(0,2*pi)),''b'')',
   'hold on',
   'plot(z,''ro'',''markerfacecolor'',''r'',''markersize'',5)',
   'axis equal',
   'axis(1.07*[-1 1 -1 1])',
   'title prevertices',
   'hold off' };
  slide(3).text={
   'There are six prevertices, but two pair are quite close together.',
   '>> z = prevertex(f); plot(z,''ro'')',
   '>> angle( z([4 1])./z([3 6] )',
   'ans =',
   '  4.8136e-006',
   '  1.1612e-005',
   '',
   'This is  the *crowding phenonenon*.  The prevertices are separated by an amount which is exponential in the effective aspect ratio of the region.'};

  %========== Slide 4 ==========

  slide(4).code={
   '' };
  slide(4).text={
   'The crowding phenomenon not only makes solution of the parameter problem difficult, but if a region locally has an aspect ratio larger than, say, 20, the prevertices may not even be distinguishable in double precision (MATLAB''s default).'};

  %========== Slide 5 ==========

  slide(5).code={
   'cla',
   'f=rectmap(p,[1 3 4 6],scmapopt(''tol'',1e-6));',
   'plot(f,4,8)',
   'title(''Map from rectangle'')' };
  slide(5).text={
   'An elegant solution to crowding is to use a canonical domain that is also elongated.  Here, we compute the map to a rectangle.  We have to specify which polygon vertices map to rectangle corners.',
   '>> f = rectmap(p,[1 3 4 6]);',
   '>> plot(f,4,8)',
   '',
   ''};

  %========== Slide 6 ==========

  slide(6).code={
   '' };
  slide(6).text={
   'The aspect ratio of the rectangle is known as the *conformal modulus*. It is determined by the geometry, and would be different for different choices of the rectangle corner images.',
   '',
   '>> modulus(f)',
   'ans =',
   '    8.8381',
   ''};

  %========== Slide 7 ==========

  slide(7).code={
   'f=stripmap(p,[3 6],scmapopt(''tol'',1e-6));',
   'plot(f,8,4)',
   'title(''Map from strip'')' };
  slide(7).text={
   'Another elongated canonical domain is the strip {0 < Im z< 1}.  This is especially useful for maps to infinite flow channels, but we can apply it here, too.',
   '>> f = stripmap(p,[3 6]);',
   '>> plot(f,8,4)',
   '',
   'Here we chose the vertices that map to the "ends" of the strip at plus/minus infinity.'};

  %========== Slide 8 ==========

  slide(8).code={
   'p = polygon([  -0.5000 - 2.5000i',
   '   0.5000 - 2.5000i',
   '   0.5000 + 1.0000i',
   '   3.0000 + 2.5000i',
   '   2.5000 + 3.0000i',
   '        0 + 1.5000i',
   '  -1.5000 + 3.5000i',
   '  -2.0000 + 3.0000i',
   '  -0.5000 + 1.0000i',
   ' ]);',
   'plot(p)' };
  slide(8).text={
   'Both the rectangle and the strip have just one primary elongation.  For regions which are multiply elongated, such as this Y, crowding again occurs.'};

  %========== Slide 9 ==========

  slide(9).code={
   'f=rectmap(p,[2 4 5 1],scmapopt(''tol'',1e-4));',
   'plot(f,4,6)' };
  slide(9).text={
   '>> f = rectmap(p,[2 4 5 1]);',
   '>> min(abs(diff(prevertex(f))))',
   '',
   'ans =',
   '    3.1196e-006',
   '',
   'No matter how one chooses two arms to match the elongation of the rectangle, the remaining arm causes crowding.'};

  %========== Slide 10 ==========

  slide(10).code={
   'f=crdiskmap(p,scmapopt(''tol'',1e-4));',
   'f=center(f,i);',
   'plot(f)',
   '' };
  slide(10).text={
   'An algorithm known as CRDT uses a representation of the prevertices that is not affected by crowding.  Using it, one can map to arbitrarily elongated regions from the disk.',
   '>> f = crdiskmap(p);',
   '>> f = center(f,i);',
   '>> plot(f)',
   '',
   'As you can see, much of the polygon maps to a tiny portion of the disk.  This is the hallmark of crowding.'};

  %========== Slide 11 ==========

  slide(11).code={
   '' };
  slide(11).text={
   'As the disk may not be a useful canonical domain for many problems, an alternative is offered. We can "rectify" the polygon by making each of its angles be an integer multiple of pi/2. Such a rectified polygon has each side parallel to one of the cartesian axes.',
   '',
   'This is accomplished by using the disk as an intermediate domain and constructing two S--C maps: the usual one, and another one with the same prevertices and rectified angles.',
   ''};

  %========== Slide 12 ==========

  slide(12).code={
   'alpha =   [ 0.5000',
   '    0.5000',
   '    1.0000',
   '    1.0000',
   '    1.5000',
   '    1.0000',
   '    1.0000',
   '    0.5000',
   '    0.5000',
   '    1.0000',
   '    1.0000',
   '    1.0000',
   '    0.5000',
   '    0.5000',
   '    1.5000',
   '    1.0000',
   '    1.0000];',
   'f = crrectmap(f,alpha);',
   'plot(rectpoly(f))',
   '' };
  slide(12).text={
   'Here we map the previous polygon to a "T" (it''s suited to that).  The angles in alpha were predetermined.',
   '>> f = crrectmap(f,alpha);',
   '',
   'Here''s a plot of the target T.',
   '>> plot(rectpoly(f))',
   '',
   'The side lengths are not freely chosen; rather, they result from the prevertices of the disk map.  This is a generalization of conformal modulus.'};

  %========== Slide 13 ==========

  slide(13).code={
   'plot(f,18,12)' };
  slide(13).text={
   'We now plot the images of vertical and horizontal lines inside the T.',
   '>> plot(f,18,12)'};
end