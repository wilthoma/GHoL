function v=readvect(cfile)
% reads a sparse vector from file produced by makemat
l=dlmread(cfile);
I= find(l(:,2)>0);
v=sparse(l(I,1)+1,l(I,2),l(I,3))';