function dX = model6(t,x,am,bm,gs,as,bs,ap,bp,k1,k2)

m = x(1); 
s = x(2); 
p = x(3);

dX1 = am - bm*m - gs*m*s ; 
dX2 = as - bs*s ; 
dX3 = ap*m*((1+k1*k2*s)/(1+k2*s)) - bp*p; 

dX = [dX1; dX2; dX3]; 