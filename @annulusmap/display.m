function display(map)
%DISPLAY   Pretty-print D-SC parameters.
%   Called automatically when a map is displayed at the command line.

%   Copyright by Alfa Heryudono, 2003.

fprintf('\n%s = DSC annulusmap:\n\n',inputname(1))
if imag(map.c) < 0
  s = '-';
else
  s = '+';
end
fprintf('  Conformal modulus = %5.10f        c = %.8g %c %.8gi \n\n', 1/map.u, real(map.c),s,abs(imag(map.c)));

lab{1} = 'Outer';
lab{2} = 'Inner';

for k=2:-1:1;
z1 = eval(['map.Z' int2str(k-1) '(:)']);
alfa1 = eval(['map.ALFA' int2str(k-1) '(:)']);
w1 = eval(['map.w' int2str(k-1) '(:)']);
phi1 = eval(['map.phi' int2str(k-1) '(:)']);
n = length(z1);

% We make disp do the heavy lifting. This way the FORMAT command works
% here too.

z1str = evalc( 'disp(z1)' );
alfa1str = evalc( 'disp(alfa1)' );
w1str = evalc( 'disp(w1)' );
phi1str = evalc( 'disp(phi1)' );

% Parse into one cell per line.
for j=1:n
  [tmp,z1str] = strtok(z1str,sprintf('\n'));  z1c{j}=tmp;
  [tmp,alfa1str] = strtok(alfa1str,sprintf('\n'));  alfa1c{j}=tmp;
  [tmp,w1str] = strtok(w1str,sprintf('\n'));  w1c{j}=tmp;
  [tmp,phi1str] = strtok(phi1str,sprintf('\n'));  phi1c{j}=tmp;
end

% Now into matrices.
z1m = strvcat(z1c);  alfa1m = strvcat(alfa1c); w1m = strvcat(w1c);  phi1m = strvcat(phi1c);

% Remove leading and trailing space blocs.
idx = find( ~all(z1m==' ') );
z1m = z1m(:,min(idx):max(idx));
idx = find( ~all(alfa1m==' ') );
alfa1m = alfa1m(:,min(idx):max(idx));
idx = find( ~all(w1m==' ') );
w1m = w1m(:,min(idx):max(idx));
idx = find( ~all(phi1m==' ') );
phi1m = phi1m(:,min(idx):max(idx));

z1v = max(size(z1m,2),6);
alfa1a = max(size(alfa1m,2),11);
w1v = max(size(w1m,2),9);
phi1a = max(size(phi1m,2),8);

b1 = blanks(2+floor((z1v-6)/2));
b2 = blanks(ceil((z1v-6)/2)+4+floor((alfa1a-11)/2));
b3 = blanks(2+floor((w1v-6)/2));
b4 = blanks(ceil((w1v-6)/2)+3+floor((phi1a-8)/2));
b5 = blanks(floor((z1v+4+alfa1a+w1v+phi1a)/2));
fprintf([b5 lab{k} ' Polygon' b5 '\n']);
fprintf( [b1 'Vertex' b2 lab{k + (-1)^(k-1)} ' angle/pi' b3 'Prevertex' b4 'Angle/pi\n'] );
%fprintf( [b1 '------' b2 '--------\n'] );

uz1v = min(size(z1m,2),6);
ualfa1a = min(size(alfa1m,2),11);
uw1v = min(size(w1m,2),6);
uphi1a = min(size(phi1m,2),11);

b1 = blanks(2+floor((6-uz1v)/2));
b2 = blanks(ceil((6-uz1v)/2)+4+floor((11-ualfa1a)/2));
b3 = blanks(floor((20-uw1v)/2));
b4 = blanks(ceil((20-uw1v)/2)+floor((5-uphi1a)/2));
str = [ repmat(b1,n,1) z1m repmat(b2,n,1) alfa1m repmat(b3,n,1) w1m repmat(b4,n,1) phi1m];

fprintf(['  ' repmat('-',1,z1v+4+alfa1a+7+w1v+4+phi1a) '\n']);
disp(str)
fprintf('\n\n')
z1m = [];alfa1m = [];w1m = [];phi1m=[];
z1c = [];alfa1c = [];w1c = [];phi1c=[];
end
