function F = rpfun(y,fdat)
%   Returns residual for solution of nonlinear equations. 
%   $Id: rpfun.m 251 2003-03-07 16:34:11Z driscoll $

[n,beta,nmlen,left,right,cmplx,qdat,corners] = deal(fdat{:});

% Transform y (unconstr. vars) to z (actual params)
z = rptrnsfm(y,corners);

% Compute the integrals appearing in nonlinear eqns.
zleft = z(left);
zright = z(right);

% To use stquad, must put strip ends into z, beta
ends = find(diff(imag(z([1:n 1]))));
z = [z(1:ends(1));Inf;z(ends(1)+1:ends(2));-Inf;z(ends(2)+1:n)];
beta = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
% Put dummy columns into qdat at ends
idx = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1];
qdat = qdat(:,[idx idx+n+1]);
% Change singularity indices to reflect ends
left = left + (left > ends(1)) + (left > ends(2));
right = right + (right > ends(1)) + (right > ends(2));

ints = zeros(size(zleft));
% Two-stage integrations
s2 = (right(:) - left(:) == 1) & (imag(zleft) - imag(zright) == 0);
mid = mean([zleft(s2) zright(s2)],2);
ints(s2) = stquadh(zleft(s2),mid,left(s2),z,beta,qdat) ...
    - stquadh(zright(s2),mid,right(s2),z,beta,qdat);

% Three-stage integrations
mid1 = real(zleft(~s2)) + i/2;
mid2 = real(zright(~s2)) + i/2;
ints(~s2) = stquad(zleft(~s2),mid1,left(~s2),z,beta,qdat) ...
    + stquadh(mid1,mid2,zeros(size(mid1)),z,beta,qdat) ...
    - stquad(zright(~s2),mid2,right(~s2),z,beta,qdat);

if any(ints==0)|any(isnan(ints))
  % Singularities were too crowded. 
  warning('Severe crowding')
end

% Compute nonlinear equation residual values.
cmplx2 = cmplx(2:length(cmplx));
F = abs(ints(~cmplx)); 		% F(1) = abs(ints(1))
F = log( (F(2:end)/F(1)) ./ nmlen(~cmplx2) );
if any(cmplx)
  F2 = log( (ints(cmplx)/ints(1)) ./ nmlen(cmplx2) );
  F = [F;real(F2);imag(F2)];
end
  