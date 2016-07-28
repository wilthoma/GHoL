function D=readmat(cfile)
%try
  M = dlmread(cfile);
  I = find(M(:,2) > 0);
  D = sparse(M(I,1),M(I,2),M(I,3))';
%catch
%fprintf('Warning: could not read %s, assuming 0\n',cfile);
%D=[];
%end