function [NoiseG, NoiseM, NoiseP] = FunTauL_TrIn(Imean,Id)
LoadFileName ='.\Model9\ModelFitIFFL9'; 
AddOn = 'VaryIndOpen'; 

Iter = 10000; 
Tau = .1;
Tmax = 150; 
NoTau = Tmax/Tau; 


load(LoadFileName); 
Gmax = 1; 
lm = .2; 
lp = lm/1.0312e5; 
Istd = 2*Imean; 

%amo = am*bm/(bm + gs*(as/bs)); 
if Id==1
r = [am; bm; gs; as; bs; ap; bp ; k1*k2]; %[am bm gs as bs ap bp k1*k2]'
else
    r = [am; bm; 0; as; bs; ap; bp ; 0];
end
x0 = [0; 0; 0; 0]; % [m s p g] 

%%                   Running Tau Leaping Algorithm 
G = zeros(Iter,1);
M = zeros(Iter,1); 
S = zeros(Iter,1); 
P = zeros(Iter,1); 


G = x0(4)*ones(Iter,1);
M = x0(1)*ones(Iter,1); 
S = x0(2)*ones(Iter,1); 
P = x0(3)*ones(Iter,1); 
I = max(0,random('norm',Imean*ones(Iter,1),Istd*ones(Iter,1))); 

% T = [0 : Tau : (NoTau-1)*Tau]; 

for i = 2:NoTau 

    Gprod = random('poiss',Tau*(lp*I.*(Gmax*ones(Iter,1)-G)));
    Gdeg = random('poiss',Tau*(lm*G));
    
    Mprod = random('poiss',Tau*(r(1)*G)); 
    Mdeg = random('poiss',Tau*r(2)*M); 
    MSdeg = random('poiss',Tau*r(3)*M.*S); 
    
    Sprod = random('poiss',Tau*(r(4)*G));
    Sdeg = random('poiss',Tau*r(5)*S);
    
    Pprod = random('poiss',Tau*(r(6)*M.*(1 + r(8)*S)));
    Pdeg = random('poiss',Tau*r(7)*P);
    
    G = min(max(0,G+ Gprod - Gdeg),Gmax); 
    M = max(0,M + Mprod - Mdeg - MSdeg); 
    S = max(0,S + Sprod - Sdeg);  
    P = max(0,P + Pprod - Pdeg);
    
end

%%                 Calculate Statistics from Tau Leaping Data 

% We are only interested in the final state (or steady state)
meanG = mean(G(:,end)); 
stdG = std(G(:,end)); 
NoiseG = stdG/meanG; 

meanM = mean(M(:,end)); 
stdM = std(M(:,end)); 
NoiseM = stdM/meanM; 

meanP = mean(P(:,end)); 
stdP = std(P(:,end)); 
NoiseP = stdP/meanP; 


end

