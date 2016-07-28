function saveashg6(v, n,deg, cfile)
f = sprintf('gra%d_%d.g6',n,deg);
gra = dlmread(f);

n=length(v);

xs = repmat(['x'],n,1);

dlmwrite(cfile, [xs f gra(1:n)]);