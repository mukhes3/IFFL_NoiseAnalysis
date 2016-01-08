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

DoxLevels = [400 666 800 1000 2000 2500 5000 15000 1000000]; 
amRec = zeros(1,length(DoxLevels)); 

for i = 1:length(DoxLevels)
    load(strcat(int2str(DoxLevels(i)),'model4params')); 
    amRec(i) = am ; 
end

DoxLevels = [400 666 800 1000 2000 2500 5000 15000 1000000]*60; 
figure(1) 

semilogx(DoxLevels/60,amRec,'LineWidth',3);
xlabel('Dox level (in ng)','FontSize',25); 
ylabel('\alpha_m','FontSize',25);
title('\alpha_m vs Dox level plot for PV','FontSize',25); 
set(gca,'FontSize',15); 

saveas(1,'AmvsDox2.jpg'); 
% saveas(1,'AmvsDox.fig'); 

%%                Fitting am against Dox Levels 
cstm = @(a,b,x)...
    a*x ./( b + x); 
sP = [171, 6.46e3]; 
[MP,gof] = fit(DoxLevels',amRec',cstm,'StartPoint', sP); 

amPred =  zeros(1,length(DoxLevels)); 
amPred2 =  zeros(1,length(DoxLevels)); 

sP2 = [171, 6.46e3, 1];

for i = 1:length(DoxLevels)
    amPred(i) = (MP.a*DoxLevels(i))/(DoxLevels(i) + MP.b); 
end

cstm2 = @(a,b,n,x)...
    a*x.^n ./( b + x.^n); 

[MP2,gof2] = fit(DoxLevels',amRec',cstm2,'StartPoint', sP2); 

for i = 1:length(DoxLevels)
    amPred2(i) = (MP2.a*(DoxLevels(i)^MP2.n))/(DoxLevels(i)^MP2.n + MP2.b); 
end

figure(2)
semilogx(DoxLevels/60,amRec,DoxLevels/60,amPred,DoxLevels/60,amPred2,'LineWidth',3);
legend('Actual','Predicted','Predicted Nonlinear','FontSize',25,'Location','BestOutside'); 
xlabel('Dox level (in ng)','FontSize',25); 
ylabel('\alpha_M','FontSize',25);
title('\alpha_M vs Dox level plot for PV','FontSize',25); 
set(gca,'FontSize',15); 

saveas(2,'AmvsDoxPredNew2.jpg'); 
% saveas(2,'AmvsDoxPredNew.fig');

save('AmPredictionNew2','amPred','MP','MP2','gof','gof2'); 