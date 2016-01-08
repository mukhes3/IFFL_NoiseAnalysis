%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Tau Leaping Algorithm Implementation 
%                                              Sumit Mukherjee     
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear; 
close all; 

%%                        Initializing Parameters 

LoadFileName ='.\IFFL_Alex_pv_March5modFISH\model4params'; 
SaveFolder ='TauL_IFFL_FISH_pvTrial_gNoise'; 
OutFileName = 'Stats.xls'; 

Iter = 2000; 
Tau = .005;
Tmax = 150; 
NoTau = Tmax/Tau; 

mkdir(SaveFolder); 
load(LoadFileName); 
Gmax = 1; 
lp = .8; 
lm = .18; 

r = [am; bm; gs; as; bs; ap; bp]; 
x0 = [65; 100; 10; 1]; % [m s p g] 

%%                   Running Tau Leaping Algorithm 
G = zeros(Iter,NoTau);
M = zeros(Iter,NoTau); 
S = zeros(Iter,NoTau); 
P = zeros(Iter,NoTau); 

G(:,1) = x0(4)*ones(Iter,1);
M(:,1) = x0(1)*ones(Iter,1); 
S(:,1) = x0(2)*ones(Iter,1); 
P(:,1) = x0(3)*ones(Iter,1); 

T = [0 : Tau : (NoTau-1)*Tau]; 

for i = 2:NoTau 

    Gprod = random('poiss',Tau*(lp*(Gmax*ones(Iter,1)-G(:,i-1))));
    Gdeg = random('poiss',Tau*(lm*G(:,i-1)));
    
    Mprod = random('poiss',Tau*(r(1)*G(:,i-1))); 
    Mdeg = random('poiss',Tau*r(2)*M(:,i-1)); 
    MSdeg = random('poiss',Tau*r(3)*M(:,i-1).*S(:,i-1)); 
    
    Sprod = random('poiss',Tau*(r(4)*G(:,i-1)));
    Sdeg = random('poiss',Tau*r(5)*S(:,i-1));
    
    Pprod = random('poiss',Tau*(r(6)*M(:,i-1)));
    Pdeg = random('poiss',Tau*r(7)*P(:,i-1));
    
    G(:,i) = min(max(0,G(:,i-1)+ Gprod - Gdeg),Gmax); 
    M(:,i) = max(0,M(:,i-1) + Mprod - Mdeg - MSdeg); 
    S(:,i) = max(0,S(:,i-1) + Sprod - Sdeg);  
    P(:,i) = max(0,P(:,i-1) + Pprod - Pdeg);
    
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

figure(1) 
plot(M'); 
xlabel('time (in hrs)'); 
ylabel('Number of molecules'); 
title('mRNA no. vs time'); 
saveas(1,strcat('.\',SaveFolder,'\M.jpg')); 

figure(2) 
plot(S'); 
xlabel('time (in hrs)'); 
ylabel('Number of molecules'); 
title('miRNA no. vs time'); 
saveas(2,strcat('.\',SaveFolder,'\S.jpg')); 

figure(3) 
plot(P'); 
xlabel('time (in hrs)'); 
ylabel('Number of molecules'); 
title('Protein no. vs time'); 
saveas(3,strcat('.\',SaveFolder,'\P.jpg')); 

figure(4) 
plot(G'); 
xlabel('time (in hrs)'); 
ylabel('Number of gene'); 
title('Activated gene no. vs time'); 
saveas(4,strcat('.\',SaveFolder,'\G.jpg')); 

save(strcat('.\',SaveFolder,'\StatsMat')); 