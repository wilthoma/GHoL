%function bestx=mymanga(A,b, x0)

alpha=.1;

mI = length(I);
% assumes m is the size and Aeq, beq havbe been set
% cou is the count
epsi=.1;
if (cou==0) 
    w=ones(mI,1);
else
    w= alpha * exp(-alpha * abs(xsI));
end
cou=cou+1
xxI=linprog([w w],[],[],Aeq(:,[I;I+m]),beq,zeros(2*mI,1),[],[],[]);
xsI=xxI(1:mI)-xxI(mI+1:2*mI);
sum(abs(xsI)>.01)






















return;

[m, n]=size(A);
alpha = 0.1;
doubleevery = 8;
toler = 1e-8;
maxiter = 100;
sptoler= .01;%=1e-5;

x=x0;
y=abs(x);
bestx=x;
bestsp = n;
for i=1:maxiter
    i
   c = alpha*exp(-alpha * y);
   f = [zeros(n,1); c];
   Aeq = [A zeros(m,n)];
   %opt = optimset('LargeScale', 'off','Simplex','off', 'Display','off');
   opt = optimset('Display','off');
   z = linprog(f, [eye(n), -eye(n); -eye(n), -eye(n)], zeros(2*n,1) ,Aeq,b,[],[],[],opt);

   x = z(1:n);
   ynew = z(n+1:end);
   if max(abs(y-ynew))<toler
       %break;
   end
   y=ynew;
   sparseness=sum(abs(x)>sptoler)
   if (sparseness<bestsp)
       bestx=x;
       bestsp=sparseness;
   end
   if mod(i, doubleevery)==0
       alpha = 2*alpha
   end
end
