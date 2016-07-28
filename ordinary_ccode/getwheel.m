function w=getwheel(n, wind)
[D DD]=readdifferential(n,0);

[m,N]=size(DD);

v= zeros(m,1);
v(wind) = 1;

beq = [zeros(N,1);1];

%w = [DD';v']\beq;

w = lsqr([DD';v'],beq,1e-15, 200);

%return;

Aeq=[DD' -DD';v' -v'];

options = optimset('LargeScale', 'off','Simplex','off');
%options = optimset();
x=linprog(ones(2*m,1),[],[],Aeq,beq,zeros(2*m,1),[],[],options);
w=x(1:m)-x(m+1:2*m);

sum(abs(DD'*w))
v'*w