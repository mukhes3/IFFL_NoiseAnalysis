function cost = modelFit10(x,sDat,mDat,pDat,tDat,am,bm,as,bs,ap,bp) 

gs = x(1); 
kf = x(2); 
kr = x(3); 
k = x(4); 
ms0 = x(5); 


m0 = mDat(1) - ms0;  
s0 = sDat(1) - ms0; 
p0 = pDat(1); 

x0 = [m0; s0; p0; ms0]; 

[Tout,Yout] = ode45(@(t,x) RunAll(t,x,am,bm,gs,kf,kr,k,ap,bp,as,bs),tDat,x0);   

function dX = RunAll(t,x,am,bm,gs,kf,kr,k,ap,bp,as,bs) 
m=x(1); 
s = x(2); 
p = x(3);
ms = x(4); 

dm = am - bm*m - kf*m*s + kr*ms; 
ds = as - bs*s - kf*m*s + kr*ms + gs*ms;
dp = ap*(m + k*ms) - bp*p;
dms = kf*m*s - kr*ms - gs*ms; 
dX = [dm; ds; dp; dms];
end

% maxM = max(Yout(:,1)); 
% maxP = max(Yout(:,3));
% cost = (norm(Yout(:,1)-mDat)^2/mDat(end)^2) + 5*(norm(Yout(:,3)-pDat)^2/pDat(end)^2) + ...
%     1*(((maxM - max(mDat))^2/mDat(end)^2) + 5*((maxP - max(pDat))^2/pDat(end)^2)) + ...
%     2*(((Yout(end,1) - mDat(end))^2/mDat(end)^2) + 5*((Yout(end,3) - pDat(end))^2/pDat(end)^2)); 

maxM = max(Yout(:,1));
maxS = max(Yout(:,2));
maxP = max(Yout(:,3));

Mn = Yout(:,1) + Yout(:,4); 
Sn = Yout(:,2) + Yout(:,4); 
P = Yout(:,3); 

cost = (norm(Mn - mDat)^2/mDat(end)^2) + (norm(Sn - sDat)^2/sDat(end)^2) + ...
       (norm(Pn - pDat)^2/pDat(end)^2) + ...
       1*(((maxM - max(mDat))^2/mDat(end)^2) + 5*((maxP - max(pDat))^2/pDat(end)^2) + ...
      ((maxS - max(sDat))^2/sDat(end)^2)) + ...
       2*(((Mn(end) - mDat(end))^2/mDat(end)^2) + 5*((Pn(end) - pDat(end))^2/pDat(end)^2) + ...
       ((Sn(end) - sDat(end))^2/sDat(end)^2)); 


end