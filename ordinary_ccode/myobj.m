function y = myobj(x, sv, solu)
y = sum(abs(sv*x + solu).^(1/2));
%y = sum(abs(sv*x + solu));
%y = sum(abs(sv*x + solu)>.01);
%I=[285 470 930 928 1296 1712 1713 734 1218 1469 598 1743 931 929 1297 1512 1470 1744]';


