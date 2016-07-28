function [D,DD,dims,nev1s,nev2s,nev3s]=findSpecSeqDims(n, deg, prefix)

[D, DD, dims] = readdifferentialhairy(n, deg, prefix);

nev1s=readNrNonEdgeVs(n,deg,prefix);
nev2s=readNrNonEdgeVs(n-1,deg-1,prefix);
nev3s=readNrNonEdgeVs(n+1,deg+1,prefix);

for nrh = 0:25 % nr of hairs
    nn = getSpecSeqConvergents(D{nrh+1}, DD{nrh+1}, nev1s{nrh+1}, nev2s{nrh+1}, nev3s{nrh+1});
    
    fprintf('%d hairs: dim=', nrh);
    fprintf(' %d,',nn);
    fprintf('\n');
end