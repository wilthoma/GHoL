function s=tomathematica(I)
s=int2str(I(1));
for i=2:length(I)
    s= strcat(s, ',',int2str(I(i)));    
end
