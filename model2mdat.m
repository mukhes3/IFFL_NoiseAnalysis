function cost = model2mdat(x,mDat,sDat,tDat,as,bs) 

am = x(1); 
bm = x(2); 
gs = x(3); 

m = mDat(1); 

J = 0; 

for i=1:5

T = tDat(i+1)-tDat(i);    

s = sDat(i);

m = am*T + (1-bm*T)*m - gs*m*s*T;
y = [m]; 
yAct = [mDat(i+1)];

J = J + (norm(y-yAct))^2; 
end

cost = J; 