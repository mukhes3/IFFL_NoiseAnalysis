%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             How aplha m varies with induction levels 
%                                                Sumit Mukherjee
%  12/16/2914 
% This curve fitting fits with respect to mean number of Inducer molecules.
% The 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear; 
close all;

%%                Comparing Dox Levels with amRec

DoxLevels = [666 800 1000 1500 2000 2500 5000 15000 1000000]; 
amRec = zeros(1,length(DoxLevels)); 

for i = 1:length(DoxLevels)
    load(strcat(int2str(DoxLevels(i)),'model4params')); 
    amRec(i) = am ; 
end

DoxLevels = [666 800 1000 1500 2000 2500 5000 15000 1000000]*60; 
figure(1) 

semilogx(DoxLevels/60,amRec,'LineWidth',3);
xlabel('Dox level (in ng)','FontSize',25); 
ylabel('\alpha_m','FontSize',25);
title('\alpha_m vs Dox level plot for PV','FontSize',25); 
set(gca,'FontSize',15); 

saveas(1,'AmvsDox.jpg'); 
% saveas(1,'AmvsDox.fig'); 

%%                Fitting am against Dox Levels 
cstm = @(a,b,x)...
    b*x ./( a + x); 

[MP,gof] = fit(DoxLevels',amRec',cstm); 

amPred =  zeros(1,length(DoxLevels)); 
amPred2 =  zeros(1,length(DoxLevels)); 

for i = 1:length(DoxLevels)
    amPred(i) = (MP.b*DoxLevels(i))/(DoxLevels(i) + MP.a); 
end

cstm2 = @(a,b,n,x)...
    b*x.^n ./( a + x.^n); 

[MP2,gof2] = fit(DoxLevels',amRec',cstm2); 

for i = 1:length(DoxLevels)
    amPred2(i) = (MP.b*DoxLevels(i)^MP2.n)/(DoxLevels(i)^MP2.n + MP2.a); 
end

figure(2)
semilogx(DoxLevels*60,amRec,DoxLevels*60,amPred,DoxLevels*60,amPred2,'LineWidth',3);
legend('Actual','Predicted','Predicted Nonlinear','FontSize',25,'Location','BestOutside'); 
xlabel('Dox level (in ng)','FontSize',25); 
ylabel('\alpha_m','FontSize',25);
title('\alpha_m vs Dox level plot for PV','FontSize',25); 
set(gca,'FontSize',15); 

saveas(2,'AmvsDoxPredNew.jpg'); 
% saveas(2,'AmvsDoxPredNew.fig');

save('AmPredictionNew','amPred','MP','MP2','gof','gof2'); 