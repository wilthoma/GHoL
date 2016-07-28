for i=1:length(I)
    newII=setdiff(1:length(I),i);
    xxx=sA(:,newII)\sb;
    if norm(sA(:,newII)*xxx-sb)<1e-6
        sum(abs(xxx))
        i
    end
end