function dX = model10(t,x,am,bm,gs,as,bs,ap,bp,kf,kr,k)

m=x(1); 
s = x(2); 
p = x(3);
ms = x(4); 

dm = am - bm*m - kf*m*s + kr*ms; 
ds = as - bs*s - kf*m*s + kr*ms + gs*ms;
dp = ap*(m + k*ms) - bp*p;
dms = kf*m*s - kr*ms - gs*ms; 
dX = [dm; ds; dp; dms];