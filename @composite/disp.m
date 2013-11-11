function disp(f)

N = length(f.maps);
for n = 1:N
  disp(sprintf('#%i',n))
  disp(f.maps{n})
end
