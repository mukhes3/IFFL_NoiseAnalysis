function dX = model2sdatcont(t,x,am,bm,gs,as,bs,ap,bp)

m = x(1); 
s = x(2); 
p = x(3);

dX1 = am - bm*m - gs*m*s ; 
dX2 = as - bs*s ; 
dX3 = ap*m - bp*p; 

dX = [dX1; dX2; dX3]; 