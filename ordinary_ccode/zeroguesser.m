
chooseN = 5;

bestsp = 1000;
bestxxx = [];
bestzerois = [];

for iter=1:10000

zerois = randi(size(sA,2), chooseN,1);

coeff=sparse(1:chooseN,zerois, ones(chooseN,1), chooseN, size(sA,2));

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