function display(p)
%   Pretty-print a polygon.

%   Copyright 1998-2003 by Toby Driscoll.
%   $Id: display.m 271 2003-05-08 18:11:36Z driscoll $

w = vertex(p);
alpha = angle(p);
n = length(w);
if n==0
    fprintf('\n empty polygon\n\n')
    return
end

fprintf('\n%s = polygon object:\n\n',inputname(1))

% We make disp do the heavy lifting. This way the FORMAT command works
% here too.

vstr = evalc( 'disp(w)' );
astr = evalc( 'disp(alpha)' );

% Parse into one cell per line.
for j=1:n
  [tmp,vstr] = strtok(vstr,sprintf('\n'));  vc{j}=tmp;
  [tmp,astr] = strtok(astr,sprintf('\n'));  ac{j}=tmp;
end

% Now into matrices.
vm = strvcat(vc);  am = strvcat(ac);

% Remove leading and trailing space blocs.
idx = find( ~all(vm==' ') );
vm = vm(:,min(idx):max(idx));
idx = find( ~all(am==' ') );
am = am(:,min(idx):max(idx));

wv = max(size(vm,2),6);
wa = max(size(am,2),8);
b1 = blanks(2+floor((wv-6)/2));
b2 =  blanks(ceil((wv-6)/2)+4+floor((wa-8)/2));
fprintf( [b1 'Vertex' b2 'Angle/pi\n'] );
%fprintf( [b1 '------' b2 '--------\n'] );

uv = min(size(vm,2),6);
ua = min(size(am,2),8);
b1 = blanks(2+floor((6-uv)/2));
b2 = blanks(ceil((6-uv)/2)+4+floor((8-ua)/2));
str = [ repmat(b1,n,1) vm repmat(b2,n,1) am ];

fprintf(['  ' repmat('-',1,wv+4+wa) '\n']);
disp(str)
fprintf('\n\n')
