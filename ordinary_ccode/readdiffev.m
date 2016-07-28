function [D DD] = readdiffev(n, deg)
%try
f = sprintf('graev%d_%d.txt',n,deg);
%f2 = sprintf('graev%d_%d.hg6',n,deg);
%if sum(size(dlmread(f2)))>0

M = dlmread(f);
if M(1,1)*M(1,2)==0  
    D=0;
else
  I = find(M(:,2) > 0);
  D = sparse(M(I,1),M(I,2),M(I,3))';
end
%else
%D=[];
%end

%catch
%fprintf('Warning: could not read %s, assuming 0\n',f);
%D=[];
%end


try
f = sprintf('graev%d_%d.txt',n+1,deg+1);
%f2 = sprintf('graev%d_%d.hg6',n+1,deg+1);

%if sum(size(dlmread(f2)))>0

M = dlmread(f);
if M(1,1)*M(1,2)==0  
    DD=0;
else
  I = find(M(:,2) > 0);
  DD = sparse(M(I,1),M(I,2),M(I,3))';
end

%else
%DD=[];
%end

catch
fprintf('Warning: could not read %s, assuming 0\n',f);
DD=[];
end

