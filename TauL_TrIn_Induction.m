%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Tau Leaping Transcription Inhibition with Induction
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear; 
close all; 

%%                        Initializing Parameters 
tStart = tic;

LoadFileName ='.\Model9\ModelFitIFFL9'; 
SaveFileName = '.\TauL_IFFL_TransInhNew\IFFLopen'; 
StatFileName = '.\TauL_IFFL_TransInhNew\StatsOpen.xls'; 
mkdir('.\TauL_IFFL_TransInhNew'); 
AddOn = '1mugOpenStd'; 

Iter = 1000; 
Tau = .1;
Tmax = 150; 
NoTau = Tmax/Tau; 


load(LoadFileName); 
Gmax = 1; 
lm = .2; 
lp = lm/1.0312e5; 
Imean = 1e6*60; 
Istd = Imean; 

%amo = am*bm/(bm + gs*(as/bs)); 
r = [am; bm; 0; as; bs; ap; bp ; 0]; %[am bm gs as bs ap bp k1*k2]'
x0 = [0; 0; 0; 0]; % [m s p g] 

%%                   Running Tau Leaping Algorithm 
G = zeros(Iter);
M = zeros(Iter); 
S = zeros(Iter); 
P = zeros(Iter); 


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

save(SaveFileName,'meanG','stdG','NoiseG','meanM','stdM','NoiseM','meanP','stdP','NoiseP','G','M','P'); 

figure(1)
set(gca,'FontSize',15);
hist(M,50); 
xlabel('mRNA counts'); 
ylabel('No. of instances'); 

figure(2)
set(gca,'FontSize',15);
hist(P,50); 
xlabel('Protein counts'); 
ylabel('No. of instances'); 

saveas(1,strcat('.\TauL_IFFL_TransInhNew\M_',AddOn,'.jpg')); 
saveas(2,strcat('.\TauL_IFFL_TransInhNew\P_',AddOn,'.jpg'));

tElapsed = toc(tStart); 

% figure(1) 
% plot(M'); 
% xlabel('time (in hrs)'); 
% ylabel('Number of molecules'); 
% title('mRNA no. vs time'); 
% saveas(1,strcat('.\',SaveFolder,'\M.jpg')); 
% 
% figure(2) 
% plot(S'); 
% xlabel('time (in hrs)'); 
% ylabel('Number of molecules'); 
% title('miRNA no. vs time'); 
% saveas(2,strcat('.\',SaveFolder,'\S.jpg')); 
% 
% figure(3) 
% plot(P'); 
% xlabel('time (in hrs)'); 
% ylabel('Number of molecules'); 
% title('Protein no. vs time'); 
% saveas(3,strcat('.\',SaveFolder,'\P.jpg')); 
% 
% figure(4) 
% plot(G'); 
% xlabel('time (in hrs)'); 
% ylabel('Number of gene'); 
% title('Activated gene no. vs time'); 
% saveas(4,strcat('.\',SaveFolder,'\G.jpg')); 
% 
% save(strcat('.\',SaveFolder,'\StatsMat')); 