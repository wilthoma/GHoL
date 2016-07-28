function dimension = getdim(nrverts, deg)
%[D,DD] = readdiffev(nrverts,deg);
[D,DD] = readdifferential(nrverts,deg);
%if size(D,

A=[D;DD'];
 dimension=min(size(A))-rank(A);
