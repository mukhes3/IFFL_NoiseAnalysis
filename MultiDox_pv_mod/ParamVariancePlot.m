%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             How aplha m varies with induction levels 
%                                                Sumit Mukherjee
%
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

figure(1) 

semilogx(DoxLevels,amRec,'LineWidth',3);
xlabel('Dox level (in ng)','FontSize',25); 
ylabel('\alpha_m','FontSize',25);
title('\alpha_m vs Dox level plot for PV','FontSize',25); 
set(gca,'FontSize',15); 

saveas(1,'AmvsDox.jpg'); 
% saveas(1,'AmvsDox.fig'); 

%%                Fitting am against Dox Levels 
cstm = @(a,b,x)...
    a*x ./( 1 + b*x); 

[MP,gof] = fit(DoxLevels',amRec',cstm); 

amPred =  zeros(1,length(DoxLevels)); 
amPred2 =  zeros(1,length(DoxLevels)); 

for i = 1:length(DoxLevels)
    amPred(i) = (MP.a*DoxLevels(i))/(MP.b*DoxLevels(i) + 1); 
end

cstm2 = @(a,b,n,x)...
    a*x.^n ./( 1 + b*x.^n); 

[MP2,gof2] = fit(DoxLevels',amRec',cstm2); 

for i = 1:length(DoxLevels)
    amPred2(i) = (MP2.a*DoxLevels(i)^MP2.n)/(MP2.b*DoxLevels(i)^MP2.n + 1); 
end

figure(2)
semilogx(DoxLevels,amRec,DoxLevels,amPred,DoxLevels,amPred2,'LineWidth',3);
legend('Actual','Predicted','Predicted Nonlinear','FontSize',25,'Location','BestOutside'); 
xlabel('Dox level (in ng)','FontSize',25); 
ylabel('\alpha_m','FontSize',25);
title('\alpha_m vs Dox level plot for PV','FontSize',25); 
set(gca,'FontSize',15); 

saveas(2,'AmvsDoxPred.jpg'); 
% % saveas(2,'AmvsDoxPred.fig');

save('AmPrediction','amPred','MP','MP2','gof','gof2'); 