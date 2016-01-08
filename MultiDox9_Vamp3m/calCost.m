function cost = calCost(x,Dox,aM)

a = x(1); 
b = x(2); 
n = x(3); 

amPred = a*Dox.^n./( b + Dox.^n); 

cost = 0*norm(amPred-aM)^2 + 10*norm(amPred(1:8)-aM(1:8))^2; 


