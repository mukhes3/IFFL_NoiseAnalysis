function dX = mpdat4(t,x,am,bm,gs,ap,bp,as,bs) 
m=x(1); 
s = x(2); 
p = x(3); 
dm = am - bm*m - gs*m*s;
ds = as - bs*s;
dp = ap*m - bp*p; 
dX = [dm; ds; dp];
end