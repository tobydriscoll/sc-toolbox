function F = dpfun(y,fdat)
%   Returns residual for solution of nonlinear equations.

[n,beta,nmlen,left,right,cmplx,qdat] = deal(fdat{:});

% Convert y values to z (prevertices)
cs = cumsum(cumprod([1;exp(-y)]));
theta = pi*cs(1:n-3)/cs(length(cs));
z = ones(n,1);
z(1:n-3) = exp(i*theta);
z(n-2:n-1) = [-1;-i];

%%% Check crowding.
%%if any(diff(theta)<eps) | any(isnan(theta))
%%  % Since abs(y) is large, use it as the penalty function.
%%  F = y;
%%  disp('Warning: Severe crowding')
%%  return
%%end

% Compute the integrals 
zleft = z(left);
zright = z(right);
angl = angle(zleft);
mid = exp(i*(angl + rem(angle(zright./zleft)+2*pi,2*pi)/2));
% For integrals between nonadjacent singularities, choose 0 as intermediate
% integration point.
mid(cmplx) = zeros(size(mid(cmplx)));
% If any are complex, the first one must be too.
if any(cmplx)
  cmplx(1) = 1;
end

ints = NaN*zleft;
ints(~cmplx) = dabsquad(zleft(~cmplx),mid(~cmplx),left(~cmplx),...
                        z,beta,qdat) + ...
    dabsquad(zright(~cmplx),mid(~cmplx),right(~cmplx),z,beta,qdat);
if any(cmplx)
  ints(cmplx) = dquad(zleft(cmplx),mid(cmplx),left(cmplx),z,beta,qdat) - ...
    dquad(zright(cmplx),mid(cmplx),right(cmplx),z,beta,qdat);
end

if any(ints==0)
  % Singularities were too crowded in practice.
%%  F = y;
  warning('Severe crowding')
end

% Compute nonlinear equation residual values.
cmplx(1) = 0;
F1 = ints(~cmplx); 		
F1 = F1(2:end)/abs(F1(1));
F2 = ints(cmplx)/ints(1);
F = [F1;real(F2);imag(F2)] - nmlen;


