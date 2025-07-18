function slide=scdfaber
% This is a slideshow file for use with playshow.m and makeshow.m
% To see it run, type 'playshow scdfaber', 

% Copyright (c) 1984-98 by The MathWorks, Inc.
% $Id: scdfaber.m 298 2009-09-15 14:36:37Z driscoll $
if nargout<1
  playshow scdfaber
else
  %========== Slide 1 ==========

  slide(1).code={
   '' };
  slide(1).text={
   ''};

  %========== Slide 2 ==========

  slide(2).code={
   '' };
  slide(2).text={
   'Let f be a conformal map from the exterior of a Jordan curve to the exterior of the unit disk such that infinity is a fixed point.  The nth Faber  polynomial is the polynomial part of the Laurent expansion of f^n at infinity.',
   '',
   'The Faber polynomials reduce to Chebyshev polynomials when the curve is a line segment, and to Taylor polynomials when it is a circle. If the curve is a polygon, the Schwarz-Christoffel exterior map can be used to compute Faber polynomial coefficients.',
   ''};

  %========== Slide 3 ==========

  slide(3).code={
   'p = polygon([  -3.0519 + 1.3846i',
   '  -0.7308 + 1.5385i',
   '  -1.8423 - 1.7538i',
   '  -1.7442 - 3.3846i',
   '   3.5000 - 1.5000i',
   '   3.0000',
   '   0.5000 - 1.0000i',
   '   2.0000 + 3.5000i',
   '  -3.0846 + 2.9231i]);',
   'plot(p)',
   'f = extermap(p);',
   'F = faber(f,20);' };
  slide(3).text={
   'We find the exterior map to this polygon.',
   '>> plot(p)',
   '>> f = extermap(p);',
   '',
   'We compute the coefficients of the first 21 Faber polynomials.',
   '>> F = faber(f,20);'};

  %========== Slide 4 ==========

  slide(4).code={
   'lim = axis;',
   '[X,Y] = meshgrid(linspace(lim(1),lim(2),50),linspace(lim(3),lim(4),50));' };
  slide(4).text={
   'Because the Faber polynomials approximate a function having unit modulus on the polygon, the lemniscate',
   '',
   '        {z: | P (z) | = 1}',
   '',
   ' for a Faber polynomial P will approximate the polygon.',
   '',
   '>> lim = axis;',
   '>> [X,Y] = meshgrid(linspace(lim(1),lim(2),50),linspace(lim(3),lim(4),50));'};

  %========== Slide 5 ==========

  slide(5).code={
   'hold on',
   'u = abs(polyval(F{6},X+i*Y));',
   '[ct,h] = contour(X,Y,u,[1 1],''r'');',
   'title(''| P_5(z) | = 1'')' };
  slide(5).text={
   '>> u = abs(polyval(F{6}, X+i*Y));',
   '>> contour(X,Y,u,[1 1],''r'')',
   '',
   'Note how badly the reentrant corners (with respect to the interior of the polygon) are represented.  In the coming frames, we will increase the degree of the Faber polynomial.'};

  %========== Slide 6 ==========

  slide(6).code={
   'delete(h)',
   'u = abs(polyval(F{11}, X+i*Y));',
   '[ct,h] = contour(X,Y,u,[1 1],''r'');',
   'title(''| P_{10}(z) | = 1'')' };
  slide(6).text={
   '>> u = abs(polyval(F{11}, X+i*Y));',
   '>> contour(X,Y,u,[1 1],''r'')'};

  %========== Slide 7 ==========

  slide(7).code={
   'delete(h)',
   'u = abs(polyval(F{16}, X+i*Y));',
   '[ct,h] = contour(X,Y,u,[1 1],''r'');',
   'title(''| P_{15}(z) | = 1'')' };
  slide(7).text={
   '>> u = abs(polyval(F{16}, X+i*Y));',
   '>> contour(X,Y,u,[1 1],''r'')'};

  %========== Slide 8 ==========

  slide(8).code={
   'delete(h)',
   'u = abs(polyval(F{21}, X+i*Y));',
   '[ct,h] = contour(X,Y,u,[1 1],''r'');',
   'title(''| P_{20}(z) | = 1'')' };
  slide(8).text={
   '>> u = abs(polyval(F{21}, X+i*Y));',
   '>> contour(X,Y,u,[1 1],''r'')',
   '',
   'Even this far, the reentrant corners are not well approximated.'};

  %========== Slide 9 ==========

  slide(9).code={
   'zp = f( exp(1i*linspace(0,2*pi)) );',
   'delete(h)',
   'title('''')' };
  slide(9).text={
   'A familiar property of Chebyshev polynomials on the interval is equioscillation.  The Faber polynomials have similar tendencies on the boundary.',
   '',
   'We create a vector of Fejer points around the polygon.',
   '>> zp = f( exp(1i*linspace(0,2*pi)) );'};

  %========== Slide 10 ==========

  slide(10).code={
   'hold off',
   'plot(1:100,abs(polyval(F{6},zp)))',
   'xlabel(''Fejer point'')',
   'ylabel(''| P_5(z) |'')' };
  slide(10).text={
   '>> plot(1:100,abs(polyval(F{6},zp)))',
   '',
   'Note the oscillations about 1.'};

  %========== Slide 11 ==========

  slide(11).code={
   'plot(1:100,abs(polyval(F{11},zp)))',
   'xlabel(''Fejer point'')',
   'ylabel(''| P_{10}(z) |'')' };
  slide(11).text={
   '>> plot(1:100,abs(polyval(F{11},zp)))',
   '',
   'Perfect equioscillation is not achieved in general, but this is fairly close.',
   ''};
end