function cost = multiDox(x,pDat,tDat,k,gs,bm,bs,ap,bp)

m0 = x(1); 
s0 = x(2);
am = x(3); 


p0 = pDat(1); 

as = k*am; 


x0 = [m0; s0; p0]; 

[Tout,Yout] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs),tDat,x0);   

function dX = mpdat4(t,x,am,bm,gs,ap,bp,as,bs) 
m=x(1); 
s = x(2); 
p = x(3); 
dm = am - bm*m - gs*m*s;
ds = as - bs*s;
dp = ap*m - bp*p; 
dX = [dm; ds; dp];
end



cost = (norm(Yout(:,3)-pDat)^2/pDat(end)^2); 
end
