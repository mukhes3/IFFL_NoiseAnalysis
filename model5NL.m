function dX = model5NL(t,x,am,bm,a,b,n,as,bs,ap,bp)

m = x(1); 
s = x(2); 
p = x(3);

dX1 = am - bm*m - (a*(m*s)^n/(1 + b*(m*s)^n)) ; 
dX2 = as - bs*s ; 
dX3 = ap*m - bp*p; 

dX = [dX1; dX2; dX3]; 