function cost = model4mpdat(x,sDat,mDat,pDat,tDat,as,bs) 

am = x(1); 
bm = x(2);
gs = x(3); 
ap = x(4); 
bp = x(5); 


m0 = mDat(1);  
s0 = sDat(1); 
p0 = pDat(1); 

x0 = [m0; s0; p0]; 

[Tout,Yout] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs),tDat,m0);   

function dX = mpdat4(t,x,am,bm,gs,ap,bp,as,bs) 
m=x(1); 
s = x(2); 
p = x(3); 
dm = am - bm*x - gs*m*s;
ds = as - bs*s;
dp = ap*m - bm*p; 
end



cost = (norm(Yout(:,1)-mDat)^2/mDat(end)^2) + (norm(Yout(:,3)-pDat)^2/pDat(end)^2); 
end
 