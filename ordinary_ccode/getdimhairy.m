function getdimhairy(nrv, deg,prefix)

if nargin<3
    prefix='';
end

[D, DD, dims]=  readdifferentialhairy(nrv, deg,prefix);

for nrh = 0:25 % nr of hairs
    A=[D{nrh+1}; DD{nrh+1}'];
    
    fprintf('%d hairs: dim(H)=%d\n', nrh, dims(nrh+1)-rank(full(A)))       
    
end