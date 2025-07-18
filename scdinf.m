function slide=scdinf
% This is a slideshow file for use with playshow.m and makeshow.m
% To see it run, type 'playshow scdinf', 

% Copyright (c) 1984-98 by The MathWorks, Inc.
% $Id: scdinf.m 298 2009-09-15 14:36:37Z driscoll $
if nargout<1
  playshow scdinf
else
  %========== Slide 1 ==========

  slide(1).code={
   '',
   '' };
  slide(1).text={
   'Press "Start" to see an explanation of infinite vertices.'};

  %========== Slide 2 ==========

  slide(2).code={
   '' };
  slide(2).text={
   ' At each vertex, a polygon has an interior angle of pi*alpha. For a finite vertex, 0 < alpha <= 2, where alpha=2 is a "crack."',
   '',
   'At an infinite vertex, there is still an interior angle in the extended complex plane (Riemann sphere).  A consistent definition is the signed angle that is swept out from the outgoing edge at the vertex, through the interior, to the incoming edge. (You can get the same result by extending the sides backward from infinity to their finite intersection.)  As a result, an infinite vertex is characterized by -2 <= alpha <= 0.',
   ''};

  %========== Slide 3 ==========

  slide(3).code={
   'axis square',
   'axis([-2 2 -2 2])',
   'box on',
   'hold on',
   'int = fill([2,0,0,2],[1,1,-1,-1],[0,.5,0],''edgecolor'',''none'');',
   'edges(1) = plot([2+i i],''color'',''b'',''linew'',2.5);',
   'edges(2) = plot([i -i],''k'',''linew'',2.5);',
   'edges(3) = plot([-i 2-i],''r'',''linew'',2.5);',
   'plot([-i,i],''o'',''markerfacecolor'',''k'',''markeredgecolor'',''k'');',
   't1 = text(-.1,1,''1/2'',''hor'',''right'',''ver'',''mid'');',
   't2 = text(-.1,-1,''1/2'',''hor'',''right'',''ver'',''mid'');',
   't3 = text(2.1,0,''\alpha=0'',''hor'',''left'',''ver'',''mid'');',
   'tl1 = text(1,-1.25,''out'',''color'',''r'',''hor'',''c'',''ver'',''t'');',
   'tl2 = text(1,1.25,''in'',''color'',''b'',''hor'',''c'',''ver'',''bot'');',
   'title(''no turn'')',
   '' };
  slide(3).text={
   'At one extreme, when alpha=0, the two edges do not subtend any interior angle at infinity.  In the plane, the adjacent edges are parallel.',
   '',
   'In this and future slides, remember that the red edge goes out to infinity, and the blue edge returns from there.  At all times the region is locally to the left of the edge as one walks along it.'};

  %========== Slide 4 ==========

  slide(4).code={
   'set(t1,''string'',''1'');',
   'set(t3,''string'',''\alpha=-1/2'');',
   'delete(tl1)',
   'delete(tl2)',
   'set(int,''xdata'',[0,0,2,2],''ydata'',[2,-1,-1,2])',
   'set(edges(1),''xd'',[0 0],''yd'',[2 1])',
   'set(edges(2),''xd'',[0 0],''yd'',[1 -1])',
   'set(edges(3),''xd'',[0 2],''yd'',[-1 -1])',
   'title(''quarter turn'')',
   '' };
  slide(4).text={
   'Here the sides take a quarter-turn when viewed at Riemannian infinity from above.'};

  %========== Slide 5 ==========

  slide(5).code={
   'set(t1,''pos'',[-0.1 0.9],''string'',''3/2'',''ver'',''top'')',
   'set(t3,''string'',''\alpha=-1'');',
   'set(int,''xdata'',[-2,0,0,2,2,-2],''ydata'',[1,1,-1,-1,2,2])',
   'set(edges(1),''xd'',[-2 0],''yd'',[1 1])',
   'set(edges(2),''xd'',[0 0],''yd'',[1 -1])',
   'set(edges(3),''xd'',[0 2],''yd'',[-1 -1])',
   'title(''half turn'')',
   '' };
  slide(5).text={
   'Now the adjacent sides are again parallel, although not necessarily collinear.  At infinity, the two edges come in from opposite directions.'};

  %========== Slide 6 ==========

  slide(6).code={
   'set(t1,''pos'',[-0.1 1.1],''string'',''7/4'',''ver'',''bot'')',
   'set(t3,''string'',''\alpha=-5/4'');',
   'set(int,''xdata'',[-2,0,0,2,2,-2],''ydata'',[-2,1,-1,-1,2,2])',
   'set(edges(1),''xd'',[-2 0],''yd'',[-2 1])',
   'set(edges(2),''xd'',[0 0],''yd'',[1 -1])',
   'set(edges(3),''xd'',[0 2],''yd'',[-1 -1])',
   'title(''five-quarters turn'')',
   '' };
  slide(6).text={
   'As alpha continues to get smaller, the region covers more of the plane.'};

  %========== Slide 7 ==========

  slide(7).code={
   'set(t1,''pos'',[0.1 0.9],''string'',''3/2'',''hor'',''left'',''ver'',''top'')',
   'set(t2,''pos'',[0.1 -0.9],''string'',''3/2'',''hor'',''left'',''ver'',''bot'')',
   'set(t3,''string'',''\alpha=-2'');',
   'set(int,''xdata'',[2,0,0,2,2,-2,-2,2],''ydata'',[-1,-1,1,1,2,2,-2,-2])',
   'set(edges(1),''xd'',[2 0],''yd'',[-1 -1])',
   'set(edges(2),''xd'',[0 0],''yd'',[1 -1])',
   'set(edges(3),''xd'',[0 2],''yd'',[1 1])',
   'title(''full turn'')',
   '' };
  slide(7).text={
   'At the final extreme, the edges are again parallel.  Now, however, they subtend a full 2pi when they meet at infinity.  We had to swap the blue and red edges in order to keep the correct sense of traveling around the boundary.'};
end