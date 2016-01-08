%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Function for Generalized Gillespie Simulation 
%                                -------   Sumit Mukherjee 
%   
%   Inputs : 
%  --------------
% 
%   Si = Input Stoichiometric Matrix 
%   So = Output Stoichiometric Matrix 
%   k  = Reaction rate matrix 
%   x0 = Initial number of molecules of each species 
%   Tmax = Final time
% 
%   Outputs : 
%  ---------------
% 
%   Xrec = Number of each specie recorded over time 
%   trec = Reaction instances 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [Xrec, trec] = GillespieGeneral(Si, So, k, x0, Tmax) 


X = x0; 
[ai,bi] = size(Si);
[ao,bo] = size(So);
Xrec1 = X;
t=0;
trec1 = t;

for i=1:5e20   %just giving a random large number so that it runs for most final times
    
    if t>=Tmax
        break;
    end
    
    a0 = 0;
    cpoints = zeros(ai,1);
    ar = zeros(ai,1);
    
    for j = 1:ai
    
    ar(j) = k(j)*prod(X.^(Si(j,:)')); 
    a0 = a0+ ar(j) ;
    cpoints(j) = sum(ar);
    
    end
    
  
    cpoints = cpoints/a0;
    
    dt=random('exp',1/a0);
        
    t=t+dt;
    
%=========Updating Concentrations of [A] and [B] after each reaction======
    
 val = random('unif',0,1);
 

temp = (cpoints>val);
[garb, rno] = max(temp); 

X = X + (So(rno,:)'-Si(rno,:)');
    
    trec1=[trec1,t];
    Xrec1 = [Xrec1 , X];
end
%=========================================================================
trec=trec1;
Xrec = Xrec1; 