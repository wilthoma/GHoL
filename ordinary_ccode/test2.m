
for iter=1:600
    
zerois = union(bestzerois,randi(length(I)));
chooseN=length(zerois);

coeff=sparse(1:chooseN,zerois, ones(chooseN,1), chooseN, length(I));

xxx=[sA; coeff] \ [sb; zeros(chooseN,1)];

if (norm(sA*xxx - sb)<1e-10 && sum(abs(xxx)) < 1e5)
    cursp = sum(abs(xxx)>.01);
    if (cursp<bestsp)
        bestsp = cursp
        bestxxx = xxx;
        bestzerois=zerois;
    end
end

end

bestsp