function cost = model2sdat(x,sDat,tDat) 

as = x(1); 
bs = x(2); 
 
s = sDat(1); 

J = 0; 

for i=1:5

T = tDat(i+1)-tDat(i);    


s = as*T + (1-bs*T)*s ;

y = [s]; 
yAct = [sDat(i+1)];

J = J + (norm(y-yAct))^2; 
end

cost = J; 