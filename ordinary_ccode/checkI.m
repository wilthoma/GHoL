I=dlmread('temp1.txt');
xxx=A(:,I)\b;
norm(A(:,I)*xxx-b)
sum(abs(xxx))
sum(abs(xxx)>.01)


size(I)