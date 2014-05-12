function  [X, info, perf] = sdogleg(fun,par, x0, opts, B0)
%SDogLeg  Secant version of Dog Leg method for nonlinear system of equations
%    f_i(x) = 0 , i=1,...,n
%  where  x  is a vector,  x = [x_1, ..., x_n] . 
%  In the discussion we also introduce the function
%    F(x) = .5 * sum(f_i(x)^2) .
%  The functions  f_i(x)  must be given by a MATLAB
%  function with declaration
%            function  f = fun(x, par)
%  par  may be dummy.
%  
%  Call:
%       [X, info {, perf}] = SDogLeg(fun,par, x0, opts {, B0})
%
%  Input parameters
%  fun  :  String with the name of the function.
%  par  :  Parameters of the function.  May be empty.
%  x0   :  Starting guess for  x .
%  opts :  Vector with six elements:
%          opts(1) = Initial trust region radius.
%          opts(2:5) used in stopping criteria:
%              ||F'||inf <= opts(2)                     or 
%              ||dx||2 <= opts(3)*(opts(3) + ||x||2)    or
%              ||f||inf <= opts(4)                      or
%              no. of iteration steps exceeds  opts(5) .
%          opts(6) = Step used in difference approximation to the 
%              Jacobian.  Not used if  B0  is present.
%
%  Output parameters
%  X    :  If  perf  is present, then array, holding the iterates
%          columnwise.  Otherwise, computed solution vector.
%  info :  Performance information, vector with 8 elements:
%          info(1:4) = final values of 
%              [||f(x)||inf  ||F'||inf  ||dx||2  Delta] 
%          info(5) = no. of iteration steps
%          info(6) = 1 :  Stopped by small  ||f(x)||inf
%                    2 :  Stopped by small  ||F'(x)||inf
%                    3 :  Stopped by small x-step
%                    4 :  Stopped by  kmax
%                    5 :  Problems, indicated by printout. 
%          info(7) = no. of function evaluations.
%          info(8) = no. of difference approximations to the Jacobian.   
%  perf :  (optional). If present, then array, holding 
%            perf(1,:) = values of  ||f(x)||inf
%            perf(2,:) = values of  ||F'(x)||inf
%            perf(3,:) = Radius of trust region,  Delta
%            perf(4,:) = values of  beta

%  Hans Bruun Nielsen,  IMM, DTU.  99.06.10

   %  Check function call
   nin = nargin;   stop = 0;
   [x n f] = check(fun,par, x0, opts, nin);   fcl = 1;
   %  Initial approximate Jacobian and inverse
   if  nin > 4
     reappr = 0;  sB = size(B0);
     if      sum(sB) == 0,  B = eye(n,n);   D = eye(n);
     elseif  any(sB ~= [n n])
       error('Dimensions of B0 do not match  x')
     else,   B = B0;  D = inv(B);  end
   else
     [B D stop] = Dapprox(fun,par,x,f,opts(6));
     reappr = 1;   fcl = fcl + n;
   end  
   Delta = opts(1);   kmax = opts(5);
   thrf = opts(4);  thrg = opts(2);  thrx = opts(3);
   nf = norm(f,inf);   F = (f'*f)/2;
   Trace = nargout > 2;
   if  Trace,  X = zeros(n,kmax+1);  perf = zeros(4,kmax+1); end 
   k = 0;   nu = 2;   nx = thrx + norm(x);
   knew = 0;   K = max(20,n);   start = 0;  beta = -1;

   while  ~stop
     %  Check simple stopping criteria
     if      nf <= thrf,  stop = 1;
     elseif  Delta <= thrx*nx,  stop = 3; 
     else
       if  reappr & (start & (nu > 16 | k-knew == K))
         % Recompute difference approximation
         [B D stop] = Dapprox(fun,par,x,f,opts(6));
         reappr = reappr + 1;   fcl = fcl + n;
         knew = k;   nu = 2;   start = 0; 
       end
     end
     if  ~stop
       %  Check gradient criterion
       g = B'*f;   ng = norm(g,inf);    k = k+1; 
       if  Trace, X(:,k) = x;   perf(:,k) = [nf ng Delta 0]'; end  
       if  ng <= thrg,  stop = 2; end
     end
     if  ~stop    %  Find step
       ng2 = norm(g);
       alpha = (ng2/norm(B*g))^2;   a = -alpha*g;   na = alpha*ng2;
       b = D*(-f);   nb = norm(b);
       if      nb <= Delta    % Newton step
         h = b;   beta = 1;   nh = nb;   dL = F;
       elseif  na >= Delta    %  Steepest descent
         h = -(Delta/ng2)*g;   beta = 0;   nh = Delta;
         dL = Delta*(ng2 - .5*Delta/alpha);
       else    % 'True' dog leg
         c = b - a;   cf = [c'*[c  2*a]  na^2-Delta^2];
         beta = max(real(roots(cf)));
         h = a + beta*c;   nh = Delta;
         dL = .5*alpha*(1-beta)^2*ng2^2 + beta*(2-beta)*F;
       end
       if  Trace,  perf(4,k) = beta; end;
       if  nh <= thrx*nx,  stop = 3;
       elseif  nh >= nx/eps
         stop = 5;   disp('Approximate Jacobian is (almost) singular')
       else
         xnew = x + h;   fn = feval(fun, xnew,par);  fcl = fcl + 1;
         Fn = (fn'*fn)/2;   dF = F - Fn;
         %  Update  B  and  D
         y = fn - f;   h = xnew - x;   nh = norm(h);
         hD = h'*D;   hDy = hD*y;
         if  hDy
           B = B + ((y - B*h)/nh^2)*h';
           D = D + ((h - D*y)/(hD*y))*hD;
         else
           stop = 5;   disp('Updating breaks down')
         end 
         if  (dL > 0) & (dF > 0)
           x = xnew;   nx = thrx + norm(x);
           F = Fn;  f = fn;   nf = norm(f,inf);   start = 1;
           Delta = Delta / max(1/3, (1 - (2*dF/dL - 1)^7));   nu = 2;
         else
           Delta = Delta / nu;  nu = 2*nu;
           if  ~start,  knew = knew+1; end
         end
       end
       if  k > kmax,  stop = 4; end 
     end
   end
   %  Set return values
   if  Trace
     X = X(:,1:k);   perf = perf(:,1:k);
   else,  X = x;  end
   info = [nf  ng  nh  Delta  k-1  stop  fcl  reappr];

% ==========  auxiliary functions  =================================

function  [x,n, f] = check(fun,par,x0,opts,nin)
%  Check function call
   sx = size(x0);   n = max(sx);
   if  (min(sx) > 1)
       error('x0  should be a vector'), end
   x = x0(:);   f = feval(fun,x,par);
   sf = size(f);
   if  any(sf ~= [n 1])
       tx = 'f  must be a column vector of the same ' 
       error([tx 'length as  x']), end
%  Thresholds
   if  nin > 4,  nopts = 5;  else,  nopts = 6; end
   if  length(opts) < nopts
       tx = sprintf('opts  must have %g elements',nopts);
       error(tx), end
   if  length(find(opts(1:nopts) <= 0))
       error('The elements in  opts  must be strictly positive'), end

function  [B,D,stop] = Dapprox(fun,par,x,f,delta)
%  Difference approximation to Jacobian and its inverse
   n = length(x);   B = zeros(n,n);
   for  j = 1 : n
     z = x;  z(j) = x(j) + delta;   d = z(j) - x(j);
     B(:,j) = (feval(fun,z,par) - f)/d;
   end 
   %  Initial  D = inv(B)  via QR factorization 
   [Q R] = qr(B);   d = abs(diag(R));
   if  min(d) < max(20,n)*eps*max(d)    
     disp('Singular approximate Jacobian')
     stop = 5;   D = [];
   else,  stop = 0;   D = inv(R) * Q'; end
       