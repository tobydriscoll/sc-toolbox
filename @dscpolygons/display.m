function display(rgn)
%   Pretty-print a doubly connected region.

%   Modification of POLYGON/DISPLAY.m
%   Modified by Alfa Heryudono, 2003.

fprintf('\n%s = doubly connected region object:\n\n',inputname(1))

ptemp = {rgn.p1,rgn.p0};

for k=1:2

p = ptemp{k};
w = vertex(p);
alpha = angle(p);
n = length(w);

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
b2 = blanks(ceil((wv-6)/2)+4+floor((wa-8)/2));

t2 = [b1 'Vertex' b2 'Angle/pi'];
bt1 = blanks(floor((length(t2)-13)/2)); % 13 is the length of thr words (Outer Polygon / Inner Polygon)
t1 = [bt1 'Inner Polygon' bt1];
uline = ['  ' repmat('-',1,wv+4+wa)];

uv = min(size(vm,2),6);
ua = min(size(am,2),8);
b1 = blanks(2+floor((6-uv)/2));
b2 = blanks(ceil((6-uv)/2)+4+floor((8-ua)/2));
str = [ repmat(b1,n,1) vm repmat(b2,n,1) am ];

if k==1
    t1temp = t1;
    t2temp = t2;
    ulinetemp = uline;
    strtemp = str;
else
    t1 = [t1temp blanks(3) [bt1 'Outer Polygon' bt1]];
    t2 = [t2temp blanks(3) t2];
    uline = [ulinetemp blanks(3) uline];
    strsize = size(str);
    if strtempsize(1) > strsize(1)
        str = [str;repmat(blanks(1),strtempsize(1)-strsize(1),strsize(2))];
        str = [strtemp repmat(blanks(3),strtempsize(1),1) str];
    elseif strtempsize(1) < strsize(1)
        strtemp = [strtemp; repmat(blanks(1),strsize(1)-strtempsize(1),strtempsize(2))];
        str = [strtemp repmat(blanks(3),n,1) str];
    else
        str = [strtemp repmat(blanks(3),n,1) str];
    end
end 
strtempsize = size(strtemp);
vm = [];am = [];
vc = [];ac = [];
end
fprintf([t1 '\n'])
fprintf([t2 '\n'])
fprintf([uline '\n'])
disp([str])
fprintf('\n\n')