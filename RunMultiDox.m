%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Generalized framework for miRNA, for multiple DOX levels
%                    Sumit Mukherjee 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear; 
close all; 

%%               Initialization 


tRange = 'C166:C174';
pRange = 'D166:D174';
datFile = 'IFFL_master_data';
OutFolderName = 'MultiDox9_Vamp3m'; 
OutFileAdd = '5000';
matFile = '.\Model9\ModelFitIFFL9'; 

tDat = xlsread(datFile,tRange);
pDat = xlsread(datFile,pRange);
load(matFile); 

pDat = pDat; 

%%                     System identification 

k = as/am; 

% Unknown parameters are [m0 s0 am] 
yMin = [0 0 0]; 
yMax = [500 300 200]; 

% k1 = 0; k2=0; 

y0 = [5 100 180]; 

options=optimset('Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolCon',1e-6,'TolFun',1e-6,'TolX',1e-6);

% y = fmincon(@(y) multiDox(y,pDat,tDat,k,gs,bm,bs,ap,bp),y0,[],[],[],[],yMin,yMax,[],options);

y = fmincon(@(y) multiDox2(y,pDat,tDat,k,gs,bm,bs,ap,bp,k1,k2),y0,[],[],[],[],yMin,yMax,[],options);

% options=gaoptimset('TolCon',1e-9,'TolFun',1e-9);
% y = ga(@(y) multiDox(y,pDat,tDat,k,gs,bm,bs,ap,bp),3,[],[],[],[],yMin,yMax,[],options);

m0 = y(1); 
s0 = y(2); 
am = y(3); 
as = am*k; 

%%                     Plotting 

TLim = 100; 

x0 = [m0; s0; pDat(1)]; 

% [Tout,Yout] = ode45(@(t,x) model2sdatcont(t,x,am,bm,gs,as,bs,ap,bp),[0 TLim],x0); 
[Tout,Yout] = ode45(@(t,x) model9(t,x,am,bm,gs,as,bs,ap,bp,k1,k2),[0 TLim],x0); 

figure(4) 
plot(Tout,Yout(:,1),'LineWidth',3);
xlabel('Time (hours)','FontSize',20); ylabel('Conc.','FontSize',20); title(strcat('Dox = ',OutFileAdd,' pg'),'FontSize',20); 
legend('m Sim','FontSize',20);

figure(5) 
plot(Tout,Yout(:,2),'LineWidth',3);
xlabel('Time (hours)','FontSize',20); ylabel('Conc.','FontSize',20); title(strcat('Dox = ',OutFileAdd,' pg'),'FontSize',20);
legend('s Sim','FontSize',20);

figure(6) 
plot(Tout,Yout(:,3),tDat,pDat,'LineWidth',3);
xlabel('Time (hours)','FontSize',20); ylabel('Conc.','FontSize',20); title(strcat('Dox = ',OutFileAdd,' pg'),'FontSize',20);
legend('p Sim','p Real','FontSize',20);

mkdir(OutFolderName); 
saveas(4,strcat('.\',OutFolderName,'\Dox',OutFileAdd,'mRNACont.jpg'));
saveas(5,strcat('.\',OutFolderName,'\Dox',OutFileAdd,'miRNACont.jpg'));
saveas(6,strcat('.\',OutFolderName,'\Dox',OutFileAdd,'pCont.jpg'));

% saveas(4,strcat('.\',OutFolderName,'\Dox',OutFileAdd,'mRNACont.fig'));
% saveas(5,strcat('.\',OutFolderName,'\Dox',OutFileAdd,'miRNACont.fig'));
% saveas(6,strcat('.\',OutFolderName,'\Dox',OutFileAdd,'pCont.fig'));

save(strcat('.\',OutFolderName,'\',OutFileAdd,'model4params'),'am','bm','gs','ap','bp','as','bs','x0','tDat','pDat');




