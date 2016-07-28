%mI = length(I);
% assumes m is the size and Aeq, beq havbe been set
% cou is the count
epsi=5.8;
mmI=size(Aeq,2)/2;
if (cou==0) 
    w=ones(1,mmI);
else
    w=1./(abs(xsI)+epsi);
end
cou=cou+1
opt = optimset('MaxIter', 500);%, 'LargeScale', 'off'); 
xxI=linprog([w w],[],[],Aeq,beq,zeros(2*mmI,1),[],[],opt);
xsI=xxI(1:mmI)-xxI(mmI+1:2*mmI);
sum(abs(xsI)>.01)
