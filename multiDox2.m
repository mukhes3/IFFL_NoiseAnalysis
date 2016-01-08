function cost = multiDox2(x,pDat,tDat,k,gs,bm,bs,ap,bp,k1,k2)

m0 = x(1); 
s0 = x(2);
am = x(3); 


p0 = pDat(1); 

as = k*am; 


x0 = [m0; s0; p0]; 

[Tout,Yout] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs,k1,k2),tDat,x0);   

function dX = mpdat4(t,x,am,bm,gs,ap,bp,as,bs,k1,k2) 
m=x(1); 
s = x(2); 
p = x(3); 
dm = am - bm*m - gs*m*s ;
ds = as - bs*s ;
dp = ap*m*(1+k1*k2*s) - bp*p; 
dX = [dm; ds; dp];
end


[pMax,pI] = max(pDat); 

cost = (norm(Yout(:,3)-pDat)^2/pI^2) + 10*(norm(Yout(end,3)-pDat(end))^2/pI^2)+...
    2*(norm(Yout(pI,3)-pMax)^2/pI^2); 
end
