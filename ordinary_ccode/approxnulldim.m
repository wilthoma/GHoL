%function r=approxnulldim(A)
[m,n]=size(A);
B=[A;sparse(randn(1,n))];
y=[zeros(m,1);1];
x=lsqr(B,y,1e-16,1000);

r=0;
solu=[];

while (norm(y-B*x)<1e-5)
   r = r+1
   solu=[solu;sparse(x')];
   %B=[sparse(1000000 * x');B];
   B=[solu;A;sparse(randn(1,n))];
   y=[zeros(m+r,1);1];
   x=lsqr(B,y,1e-16,1000);
end