

for iter=1:700
redIbak=redI;
i=randi(length(redI));
redI(i)=[];
xxx = A(:,redI)\b;
relres=norm(A(:,redI)*xxx-b);
if (relres>1e-5)
  redI=redIbak;
  fprintf('no :(\n');
else
  fprintf('found one\n');
end

end

