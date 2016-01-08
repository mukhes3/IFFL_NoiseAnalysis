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

DoxLevels = [400 666 800 1000 1500 2000 2500 5000 15000 1000000]; 
amRec = zeros(1,length(DoxLevels)); 

for i = 1:length(DoxLevels)
    load(strcat(int2str(DoxLevels(i)),'model4params')); 
    amRec(i) = am ; 
end

DoxLevels = [400 666 800 1000 1500 2000 2500 5000 15000 1000000]*60; 
figure(1) 

semilogx(DoxLevels/60,amRec,'LineWidth',3);
xlabel('Dox level (in ng)','FontSize',25); 
ylabel('\alpha_m','FontSize',25);
title('\alpha_m vs Dox level plot for PV','FontSize',25); 
set(gca,'FontSize',15); 

saveas(1,'AmvsDox2.jpg'); 
% saveas(1,'AmvsDox.fig'); 

%%                Fitting am against Dox Levels 

options=optimset('Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolCon',1e-6,'TolFun',1e-6,'TolX',1e-6);
yMin = [100 .1e3 1]; 
yMax = [300 1e4 1]; 
sP = [171 6.46e3 1]; 

[y,fval] = fmincon(@(y) calCost(y,DoxLevels,amRec),sP,[],[],[],[],yMin,yMax,[],options);

a1 = y(1); 
b1 = y(2); 

amPred =  zeros(1,length(DoxLevels)); 
amPred2 =  zeros(1,length(DoxLevels)); 

for i = 1:length(DoxLevels)
    amPred(i) = (a1*DoxLevels(i))/(DoxLevels(i) + b1); 
end

yMin = [0 0 1]; 
yMax = [300 1e4 4]; 
sP2 = [171 6.46e3 1]; 

[y2,fval2] = fmincon(@(y) calCost(y,DoxLevels,amRec),sP,[],[],[],[],yMin,yMax,[],options);

a2 = y2(1); 
b2 = y2(2); 
n2 = y2(3); 

for i = 1:length(DoxLevels)
    amPred2(i) = (a2*(DoxLevels(i)^n2))/(DoxLevels(i)^n2 + b2); 
end

figure(2)
semilogx(DoxLevels/60,amRec,DoxLevels/60,amPred,DoxLevels/60,amPred2,'LineWidth',3);
legend('Actual','Predicted','Predicted Nonlinear','FontSize',25,'Location','BestOutside'); 
xlabel('Dox level (in ng)','FontSize',25); 
ylabel('\alpha_M','FontSize',25);
title('\alpha_M vs Dox level plot for PV','FontSize',25); 
set(gca,'FontSize',15); 

saveas(2,'AmvsDoxPredNew3.jpg'); 
% saveas(2,'AmvsDoxPredNew.fig');

save('AmPredictionNew3','amPred','a1','b1','a2','b2','n2','fval','fval2'); 