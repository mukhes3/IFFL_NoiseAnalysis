function cost = model4sdat(x,sDat,tDat) 

as = x(1); 
bs = x(2); 
 
s0 = sDat(1); 

[Tout,Yout] = ode45(@(t,x) sdat4(t,x,as,bs),tDat,s0);   

function ds = sdat4(t,x,as,bs)
ds = as - bs*x; 
end



cost = norm(Yout-sDat)^2; 
end
 