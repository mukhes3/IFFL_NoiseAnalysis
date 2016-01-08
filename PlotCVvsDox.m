%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Creates plot of Dox vs CV plot 
% Id = 1 means closed loop 
% Id = 0 means open loop 
%   12/16/2014 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear; 
close all; 

%%                       Definiing parameters 

Dox = [400 666 800 1000 1500 2000 2500 5000 15000 1000000];

l = length(Dox); 

nG = zeros(1,l); 
nM = zeros(1,l);
nP = zeros(1,l);

nGo = zeros(1,l); 
nMo = zeros(1,l);
nPo = zeros(1,l);

%Closed Loop
for i = 1:l 
    [nG(i),nM(i),nP(i)] = FunTauL_TrIn(Dox(i)*60,1); 
end

%Open Loop 
for i = 1:l 
    [nGo(i),nMo(i),nPo(i)] = FunTauL_TrIn(Dox(i)*60,0); 
end
%%                     Plot Figure 

figure(1) 
set(gca,'FontSize',15);
semilogx(Dox,nP,Dox,nPo,'LineWidth',3)
xlabel('Dox'); ylabel('CV_p'); 
legend('IFFL','Open Loop'); 

figure(2) 
set(gca,'FontSize',15);
semilogx(Dox,nM,Dox,nMo,'LineWidth',3)
xlabel('Dox'); ylabel('CV_m'); 
legend('IFFL','Open Loop'); 



saveas(1,'.\ExtIntNoise\CVpVsDox.jpg'); 
saveas(2,'.\ExtIntNoise\CVmVsDox.jpg'); 