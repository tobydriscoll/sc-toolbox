function [F,M] = faber(M,m)
%FABER  Faber polynomial coefficients for polygonal regions.
%   F = FABER(M,K) computes the coefficients of the Faber polynomials up
%   to degree K for the polygon in the extermap M. F will be a cell
%   array of length K+1 such that F{j} is a vector of coefficients for
%   the Faber polynomial of degree j-1. Note that the leading (highest
%   degree) coefficient is always first, to be compatible with POLYVAL.
% 
%   [F,M] = FABER(P,K) first computes an extermap for the polygon P.
%   
%   See also EXTERMAP, DEMOFAB.

%   Copyright 1998 by Toby Driscoll.
%   $Id: faber.m 298 2009-09-15 14:36:37Z driscoll $

%   This function follows somewhat closely the procedure outlined in
%   section 4 of Starke and Varga (Num. Math., 1993), except that no
%   symmetry of the polygon is assumed.

% Parse input
if isa(M,'extermap')
  p = polygon(M);
elseif isa(M,'polygon')
  p = M;
  M = extermap(p);
else
  msg = sprintf('Expected %s to be of class extermap or polygon.',...
      inputname(1));
  error(msg)
end  

% Low-level quantities
beta = 1-angle(p);
r = parameters(M);
z = r.prevertex;
c = r.constant;

n = length(p);
gam = ones(n,m+1);			% coeffs of binomial expansion
for k = 1:m
  gam(:,k+1) = -gam(:,k).*(beta-k+1)./(k*z);
end

% Compute the coeffs of the Laurent expansion of Psi
e1 = zeros(m+1,1);
e1(1) = 1;
C = -c*e1;
for j = 1:n
  C = toeplitz(gam(j,:).', e1')*C;
end
C = C(3:m+1)./(-(1:m-1)');
x0 = 10^(-10/m);
c0 = eval(M,x0) + c/x0 - x0.^(1:m-1)*C;
C = [c0;C];

% Use the Faber recurrence to compute polynomial coeffs
P = zeros(m+1,m+1);
P(1,1) = 1;				% poly coeffs, low order first
F = cell(1,m+1);			% high order first (MATLAB style)
F{1} = 1;
for k = 1:m
  P(1:k+1,k+1) = ([0;P(1:k,k)] - P(1:k+1,1:k)*[k*C(k);C(k-1:-1:1)])/(-c);
  F{k+1} = flipud(P(1:k+1,k+1));
  % Normalize so that leading coefficients are real
  F{k+1} = F{k+1} / sign(F{k+1}(1));
end
