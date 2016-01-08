function cost = model9mp(x,sDat,mDat,pDat,tDat,as,bs) 

am = x(1); 
bm = x(2);
gs = x(3); 
k1 = x(4);
k2 = x(5); 
ap = x(6); 
bp = x(7); 


m0 = mDat(1);  
s0 = sDat(1); 
p0 = pDat(1); 

x0 = [m0; s0; p0]; 

[Tout,Yout] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,k1,k2,ap,bp,as,bs),tDat,x0);   

function dX = mpdat4(t,x,am,bm,gs,k1,k2,ap,bp,as,bs) 
m=x(1); 
s = x(2); 
p = x(3); 
dm = am*(1 + k2*s) + (k2*m*(as-bs*s)/(1 + k2*s))- bm*m - gs*m*s ; 
ds = as - bs*s;
dp = ap*m*((1+k1*k2*s)/(1+k2*s)) - bp*p; 
dX = [dm; ds; dp];
end


maxM = max(Yout(:,1)); 
maxP = max(Yout(:,3));
cost = (norm(Yout(:,1)-mDat)^2/mDat(end)^2) + 5*(norm(Yout(:,3)-pDat)^2/pDat(end)^2) + ...
    1*(((maxM - max(mDat))^2/mDat(end)^2) + 5*((maxP - max(pDat))^2/pDat(end)^2)) + ...
    2*(((Yout(end,1) - mDat(end))^2/mDat(end)^2) + 5*((Yout(end,3) - pDat(end))^2/pDat(end)^2)); 
end


