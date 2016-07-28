
sD=D2(:,I)'*D2(:,I);
sA=[sD;w57(I)';w39(I)'];
sb=[zeros(size(sD,1),1);1;0];
xxx=sA\sb;
norm(sA*xxx-sb)
sum(abs(xxx)>.01)

