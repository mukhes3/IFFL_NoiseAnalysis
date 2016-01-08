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

mRange = 'D146:D151';
sRange = 'D26:D31'; 
tRange = 'C26:C31';
pRange = 'D32:D37';
datFile = 'AlexIFFLdata.xls';
savName = 'Model9_p'; 


mDat = xlsread(datFile,mRange); 
sDat = xlsread(datFile,sRange);
tDat = xlsread(datFile,tRange);
pDat = xlsread(datFile,pRange);

%%                     Parameter Estimation  for miRNA

xMin = [40 0.025]; 
xMax = [60 0.035]; 
% xMax = 1e3*ones(4,1); 
% phi = [am bm gs as bs];

x0 = [ 50.002 0.0264]; 

options=optimset('Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolCon',1e-9,'TolFun',1e-9,'TolX',1e-9);

x = fmincon(@(x) model4sdat(x,sDat,tDat),x0,[],[],[],[],xMin,xMax,[],options); 

 
as = x(1); 
bs = x(2); 


%%                  Parameter estimation for mRNA 
%[am bm gs ap bp] 

yMin = [100 0.03 0 0 0.01]; %[200 .001 0 0 .001]
yMax = [300 0.5 0 10 0.10]; %[3000 .6 .001 5000 1]

y0 = [140 0.05 0 3.6378 0.0408]; %[248 0.19 3.2256e-2 450 .03] 

y = fmincon(@(y) model4mpdat(y,sDat,mDat,pDat,tDat,as,bs),y0,[],[],[],[],yMin,yMax,[],options); 

am = y(1); 
bm = y(2); 
gs = y(3); 
ap = y(4); 
bp = y(5); 


%  am = 248.5884; 
%  bm = .185;
%  gs = 3.0256e-4; 


%%                Running ODE model 
% ap = 1.7186; 
% bp = .03; 


TLim = 100; 

x0 = [mDat(1); sDat(1); pDat(1)]; 

[Tout,Yout] = ode45(@(t,x) model2sdatcont(t,x,am,bm,gs,as,bs,ap,bp),[0 TLim],x0); 

figure(4) 
set(gca,'FontSize',20);
plot(Tout,Yout(:,1),tDat,mDat,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('m Sim','m Real','Best');

figure(5) 
set(gca,'FontSize',20);
plot(Tout,Yout(:,2),tDat,sDat,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('s Sim','s Real','Best');

figure(6) 
set(gca,'FontSize',20);
plot(Tout,Yout(:,3),tDat,pDat,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('p Sim','p Real','Location','Best');
mkdir(savName); 
saveas(4,strcat('.\',savName,'\model4mRNACont.jpg'));
saveas(5,strcat('.\',savName,'\model4miRNACont.jpg'));
saveas(6,strcat('.\',savName,'\model4pCont.jpg'));


save(strcat('.\',savName,'\model4params'),'am','bm','gs','ap','bp','as','bs');

%% Optional Section (plotting open loop with current parameters) 







