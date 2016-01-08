function cost = model2pdat(x,mDat, pDat,tDat) 

ap = x(1); 
bp = x(2); 
 
p = pDat(1); 

J = 0; 

for i=1:5

T = tDat(i+1)-tDat(i);    


p = ap*T*mDat(i) + (1-bp*T)*p ;

y = [p]; 
yAct = [pDat(i+1)];

J = J + (norm(y-yAct))^2; 
end

cost = J; 