%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%             Model Fitting IFFL with ODE45 
%                                  Sumit Mukherjee 
%
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear; 
close all; 

%%                      Initialization

mRange = 'G2:G7';
sRange = 'D8:D13'; 
tRange = 'C2:C7';
pRange = 'D14:D19';
datFile = 'AlexIFFLdata.xls';
SavDir = 'Model10'; 
mkdir(SavDir); 

mDat = xlsread(datFile,mRange); 
sDat = xlsread(datFile,sRange);
tDat = xlsread(datFile,tRange);
pDat = xlsread(datFile,pRange);

%%                     Parameter Estimation  for All

am = 140.7183; 
ap = 6.5936; 
as = 51.0592; 
bm = .3619; 
bp = .0201; 
bs = .0247; 

%x = [ gs kf kr k ms0] 

xMin = [0 0 0 0 0 ]; 
xMax = [100 1000 1000 .99 mDat(1)]; 
% xMax = 1e3*ones(4,1); 
% phi = [am bm gs as bs];

x0 = [ 1 100 200 .5 mDat(1)*.01]; 

options=optimset('Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolCon',1e-6,'TolFun',1e-6,'TolX',1e-6);

x = fmincon(@(x) modelFit10(x,sDat,mDat,pDat,tDat,am,bm,as,bs,ap,bp),x0,[],[],[],[],xMin,xMax,[],options); 

 
gs = x(1); 
kf = x(2); 
kr = x(3); 
k = x(4); 
ms0 = x(5);


%%                Running ODE model 
% ap = 1.7186; 
% bp = .03; 


TLim = 100; 

x0 = [mDat(1)-ms0; sDat(1)-ms0; pDat(1); ms0]; 

[Tout,Yout] = ode45(@(t,x) model9(t,x,am,bm,gs,as,bs,ap,bp,k1,k2),[0 TLim],x0); 

figure(4) 
set(gca,'FontSize',15);
plot(Tout,(Yout(:,1)+Yout(:,4)),tDat,mDat,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('m Sim','m Real');

figure(5) 
set(gca,'FontSize',15);
plot(Tout,(Yout(:,2)+Yout(:,4)),tDat,sDat,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('s Sim','s Real');

figure(6) 
set(gca,'FontSize',15);
plot(Tout,Yout(:,3),tDat,pDat,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('p Sim','p Real');


saveas(4,strcat('.\',SavDir,'\model10mRNA.jpg'));
saveas(5,strcat('.\',SavDir,'\model10miRNA.jpg'));
saveas(6,strcat('.\',SavDir,'\model10protein.jpg'));

save(strcat('.\',SavDir,'\ModelFitIFFL10'),'am','bm','gs','k1','k2','ap','bp','as','bs');


