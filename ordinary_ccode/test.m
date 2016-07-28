
rrI=[];
[n,m]=size(sD);

fsD=full(sD);
%initrk=rank(fsD);

for i=1:m
    i
    v=zeros(1,m);
    v(i)=1;
    %newrk=rank([fsA;v]);
    mxx=[sD;v]\[zeros(n,1);1];
    if (norm([sD;v]*mxx-[zeros(n,1);1])<1e-8 && sum(abs(mxx))<1e5)
    
    %if newrk>initrk
        rrI=[rrI;i];
        size(rrI)
    end
end
    