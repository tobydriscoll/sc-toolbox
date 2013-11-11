function w = feval(f,z)
%FEVAL  Evaluate a composite map.

w = z;
for n = 1:length(f.maps)
  w = feval(f.maps{n},w);
end
