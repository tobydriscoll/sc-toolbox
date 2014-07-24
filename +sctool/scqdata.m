function qdat = scqdata(beta,nqpts);
%SCQDATA Gauss-Jacobi quadrature data for SC Toolbox.
%   SCQDATA(BETA,NQPTS) returns a matrix of quadrature data suitable for
%   other SC routines.  BETA is a vector of turning angles corresponding
%   to *finite* singularities (prevertices and, for exterior map, the
%   origin).  NQPTS is the number of quadrature points per subinterval,
%   roughly equal to -log10(error).
%   
%   All the SC routines call this routine as needed, and the work
%   required is small, so you probably never have to call this function
%   directly.
%   
%   See also GAUSSJ, HPPARAM, DPARAM, DEPARAM, STPARAM, RPARAM.

%   Copyright 1998 by Toby Driscoll.
%   $Id: scqdata.m 298 2009-09-15 14:36:37Z driscoll $

n = length(beta);
qnode = zeros(nqpts,n+1);
qwght = zeros(nqpts,n+1);
for j = find(beta(:)>-1)'
  [qnode(:,j),qwght(:,j)] = gaussj(nqpts,0,beta(j));
end
[qnode(:,n+1),qwght(:,n+1)] = gaussj(nqpts,0,0);
qdat = [qnode,qwght];
