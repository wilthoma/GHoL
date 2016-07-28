function Nrs = readNrNonEdgeVs(n,deg,prefix)

if nargin<3
    prefix = '';
end

f = sprintf('gra%shairy%d_%d.nev',prefix,n,deg);

MAX_HAIRS=30;

Nrs=cell(MAX_HAIRS+1,1);

% skip empty files
c=dir(f);
if c.bytes == 0
  for nrh=0:MAX_HAIRS 
     Nrs{nrh+1} =[]; 
  end
  return;
end
MM = dlmread(f);


for nrh=0:MAX_HAIRS
  J=find(MM(:,1)==nrh);
  if (isempty(J))
      Nrs{nrh+1} = [];    
  else
      Nrs{nrh+1} =MM(J, 2);
  end
end