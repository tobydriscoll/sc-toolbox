function slide=scdtutor
% This is a slideshow file for use with playshow.m and makeshow.m
% To see it run, type 'playshow scdemtut', 

% Copyright (c) 1984-98 by The MathWorks, Inc.
% $Id: scdtutor.m 298 2009-09-15 14:36:37Z driscoll $
if nargout<1,
  playshow scdtutor
else
  %========== Slide 1 ==========

  slide(1).code={
   '' };
  slide(1).text={
   'Press the "Start" button to see a tutorial demonstration of the elementary features of the SC Toolbox.'};

  %========== Slide 2 ==========

  slide(2).code={
   'p = polygon([i; -1+i; -1-i; 1-i; 1; 0]);',
   'plot(p)',
   '' };
  slide(2).text={
   'We begin with an old MATLAB standby, the L-shaped region.',
   '>> p = polygon([i; -1+i; -1-i; 1-i; 1; 0]);',
   '>> plot(p)',
   ''};

  %========== Slide 3 ==========

  slide(3).code={
   'opt = sctool.scmapopt(''trace'',''off'');',
   'f = hplmap(p,opt);',
   'plot(f,12,6)',
   '' };
  slide(3).text={
   'Now we solve the parameter problem for a map from the half-plane, and visualize the map by plotting the image of a square grid of 12 vertical and 6 horizontal lines.',
   '>> f = hplmap(p);',
   '>> plot(f,12,6)',
   ''};

  %========== Slide 4 ==========

  slide(4).code={
   '' };
  slide(4).text={
   'Note how the lines intersect at right angles.  Also note how they converge at the last vertex (the origin), since by default that is the image of the point at infinity.'};

  %========== Slide 5 ==========

  slide(5).code={
   'f = diskmap(f);',
   'plot(f)' };
  slide(5).text={
   'Let''s change the fundamental domain from the half-plane to the unit disk.',
   '>> f = diskmap(f);',
   '>> plot(f)'};

  %========== Slide 6 ==========

  slide(6).code={
   'f = center(f,0.4-0.6i);',
   'plot(f)',
   '' };
  slide(6).text={
   'The radial curves intersect at the point which is the image of zero, usually called the conformal center.  We can reset the center to be wherever we like.',
   '>> f = center(f,0.4-0.6i);',
   '>> plot(f)',
   ''};

  %========== Slide 7 ==========

  slide(7).code={
   'plot(exp(i*pi/4)*(f-center(f)))' };
  slide(7).text={
   'The map f and the polygon p are both MATLAB objects.  Each encapsulates all necessary data and has certain functions and operations defined for it, some of which extend the definitions of standard MATLAB functions.  That''s how we can make sense of statements such as "plot(f)".',
   '',
   'We can also apply affine transformations to a polygon or map.',
   '>> plot( exp(i*pi/4)*(f-center(f)) )'};

  %========== Slide 8 ==========

  slide(8).code={
   'plot(exp(i*linspace(0,2*pi))), axis square, axis(1.05*[-1 1 -1 1])',
   'hold on, plot(prevertex(f),''ro'')',
   'title(''Prevertices of the disk map'')',
   'hold off' };
  slide(8).text={
   'If you want to extract the data that is hidden in the object, there are functions to do that.',
   '>> plot(exp(i*linspace(0,2*pi)))',
   '>> plot(prevertex(f),''ro'')',
   ''};

  %========== Slide 9 ==========

  slide(9).code={
   'cla,drawnow',
   'f = extermap(p,opt);',
   'plot(f)' };
  slide(9).text={
   'We can also map from the unit disk to the exterior of the L.',
   '>> f = extermap(p,opt);',
   '>> plot(f)',
   '',
   'You may think of these as equipotential and flow lines caused by an L-shaped conductor.'};

  %========== Slide 10 ==========

  slide(10).code={
   'f = rectmap(p,[1 2 4 5],opt);',
   'plot(f)',
   'title(sprintf(''Resistance =  %.4f'',modulus(f)))' };
  slide(10).text={
   'Another important type of fundamental domain is the rectangle.  We have to choose which vertices will map to the corners of the rectangle.  That determines the aspect ratio of the rectangle uniquely.',
   '',
   '>> f = rectmap(p,[1 2 4 5],opt);',
   '>> plot(f)',
   '',
   'You may notice that these computations are slower than the previous cases. That''s because the SC integrand now involves hyperbolic functions.'};

  %========== Slide 11 ==========

  slide(11).code={
   '' };
  slide(11).text={
   'These are the equipotentials and flow lines in a resistor.  Because resistance is a conformal invariant, the resistance (in suitable units) is just the aspect ratio of the source rectangle.'};

  %========== Slide 12 ==========

  slide(12).code={
   'cla',
   'f = stripmap(p,[1 4],opt);',
   'plot(f)' };
  slide(12).text={
   'Still another type of fundamental domain is the infinite strip {| Im z | < 1}.  We get to choose the images of the ends of the strip.  Typically these woud be the ends of an infinite  polygonal channel,  but that is not a requirement.',
   '>> f = stripmap(p,[1 4],opt);',
   '>> plot(f)',
   '',
   'Notice how the lines converge at the images of infinity.'};

end