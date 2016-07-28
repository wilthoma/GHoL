function [D DD dim] = readdifferentialhairy(n, deg,prefix)

if nargin<3
    prefix = '';
end

%try
f = sprintf('gra%shairy%d_%d.txt',prefix,n,deg);
%f2 = sprintf('grahairy%d_%d.hg6',n,deg);
%if sum(size(dlmread(f2)))>0
MAX_HAIRS=30;

MM = dlmread(f);
D=cell(MAX_HAIRS+1,1);
dim=zeros(MAX_HAIRS+1,1);
for nrh=0:MAX_HAIRS
  J=find(MM(:,1)==nrh);  
  M=MM(J, 2:end);
  if isempty(M)
      D{nrh+1}=[];
      continue
  end
  dim(nrh+1) = max(M(:,1));
  if M(1,1)*M(1,2)==0 && size(M,1)==1 
    D{nrh+1}=[];
  else
    I = find(M(:,2) > 0);
    D{nrh+1} = sparse(M(I,1),M(I,2),M(I,3))';
  end
end
%else
%D=[];
%end

%catch
%fprintf('Warning: could not read %s, assuming 0\n',f);
%D=[];
%end


try
f = sprintf('gra%shairy%d_%d.txt',prefix,n+1,deg+1);
%f2 = sprintf('gra%d_%d.hg6',n+1,deg+1);

%if sum(size(dlmread(f2)))>0
DD=cell(MAX_HAIRS+1,1);
MM = dlmread(f);
for nrh=0:MAX_HAIRS
  J=find(MM(:,1)==nrh);
  M=MM(J, 2:end);
  if isempty(M)
      DD{nrh+1}=[];
      continue
  end
  if M(1,1)*M(1,2)==0 && size(M,1)==1 
    DD{nrh+1}=[];
  else
    I = find(M(:,2) > 0);
    DD{nrh+1} = sparse(M(I,1),M(I,2),M(I,3))';
  end
end
%else
%DD=[];
%end

catch
fprintf('Warning: could not read %s, assuming 0\n',f);
DD=[];
end

