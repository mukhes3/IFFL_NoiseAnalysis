function cost = model1(x,mDat,sDat,tDat) 

am = x(1); 
bm = x(2); 
gs = x(3); 
as = x(4); 
bs = x(5); 

m = mDat(1); 
s = sDat(1); 

J = 0; 

for i=1:5

T = tDat(i+1)-tDat(i);    

m = am*T + (1-bm*T)*m - gs*m*s*T; 
s = as*T + (1-bs*T)*s ;

y = [m;s]; 
yAct = [mDat(i+1);sDat(i+1)];

J = J + (norm(y-yAct))^2; 
end

cost = J; 