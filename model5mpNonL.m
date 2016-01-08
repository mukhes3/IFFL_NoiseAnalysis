function cost = model5mpNonL(x,sDat,mDat,pDat,tDat,as,bs) 

am = x(1); 
bm = x(2);
a = x(3);
b = x(4); 
n = x(5); 
ap = x(6); 
bp = x(7); 


m0 = mDat(1);  
s0 = sDat(1); 
p0 = pDat(1); 

x0 = [m0; s0; p0]; 

[Tout,Yout] = ode45(@(t,x) mpdat4(t,x,am,bm,a,b,n,ap,bp,as,bs),tDat,x0);   

function dX = mpdat4(t,x,am,bm,a,b,n,ap,bp,as,bs) 
m=x(1); 
s = x(2); 
p = x(3); 
dm = am - bm*m - (a*(m*s)^n/(1 + b*(m*s)^n));
ds = as - bs*s;
dp = ap*m - bp*p; 
dX = [dm; ds; dp];
end



cost = (norm(Yout(:,1)-mDat)^2/mDat(end)^2) + 5*(norm(Yout(:,3)-pDat)^2/pDat(end)^2); 
end
