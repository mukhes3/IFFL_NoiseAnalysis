function cost = model4pdat(x,mDat, pDat,tDat) 

ap = x(1); 
bp = x(2); 
 
p0 = pDat(1);  

[Tout,Yout] = ode45(@(t,x) pdat4(t,x,ap,bp,mDat,tDat),tDat,p0);   

function dp = pdat4(t,x,ap,bp,mDat,tDat) 
    
m = mDat(I); 
dp = ap*m - bp*x; 
end



cost = norm(Yout-pDat)^2; 
end
 