mI = length(I);
m = size(Aeq,2)/2;
% assumes m is the size and Aeq, beq havbe been set
% cou is the count
epsi=.3;
if (cou==0) 
    w=ones(mI,1);
else
    w=1./(abs(xsI)+epsi);
end
cou=cou+1
xxI=linprog([w w],[],[],Aeq(:,[I;I+m]),beq,zeros(2*mI,1),[],[],[]);
xsI=xxI(1:mI)-xxI(mI+1:2*mI);
sum(abs(xsI)>.01)
