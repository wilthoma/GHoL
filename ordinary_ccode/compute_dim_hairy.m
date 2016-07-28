% computes dimensions of hairy graph homology

%nrv = 10 % nr vertices
nrv=14
deg = 12

[D, DD, dims]=  readdifferentialhairy(nrv, deg);
dims
for nrh = 0:25 % nr of hairs
    A=[D{nrh+1}; DD{nrh+1}'];
    
    fprintf('%d hairs: dim(H)=%d\n', nrh, dims(nrh+1)-rank(full(A)))       
    
end

