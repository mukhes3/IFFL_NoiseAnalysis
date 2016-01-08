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

mRange = 'D158:D166';
sRange = 'D167:D175'; 
tRange = 'C167:C175';
pRange = 'D176:D184';
datFile = 'AlexIFFLdata.xls';
SavDir = 'Model9'; 
mkdir(SavDir); 

mDat = xlsread(datFile,mRange); 
sDat = xlsread(datFile,sRange);
tDat = xlsread(datFile,tRange);
pDat = xlsread(datFile,pRange);

%%                     Parameter Estimation  for miRNA

xMin = [0 .02]; 
xMax = [500 .04]; 
% xMax = 1e3*ones(4,1); 
% phi = [am bm gs as bs];

x0 = [ 45 .03]; 

options=optimset('Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolCon',1e-6,'TolFun',1e-6,'TolX',1e-6);

x = fmincon(@(x) model4sdat(x,sDat,tDat),x0,[],[],[],[],xMin,xMax,[],options); 

 
as = x(1); 
bs = x(2); 


%%                  Parameter estimation for mRNA

%[ am bm gs k1 k2 ap bp] 
yMin = [100 .25 0 1e-8 1e-5 0.1 .01]; 
yMax = [3000 .5 .01 .9999 1000 20 .1];

y0 = [140 .36 0.0018 1e-7 1e-3 6.59 .02]; 

y = fmincon(@(y) model9mp(y,sDat,mDat,pDat,tDat,as,bs),y0,[],[],[],[],yMin,yMax,[],options); 

am = y(1); 
bm = y(2); 
gs = y(3);
k1 = y(4); 
k2 = y(5); 
ap = y(6); 
bp = y(7); 


%  am = 248.5884; 
%  bm = .185;
%  gs = 3.0256e-4; 


%%                Running ODE model 
% ap = 1.7186; 
% bp = .03; 


TLim = 100; 

x0 = [mDat(1); sDat(1); pDat(1)]; 

[Tout,Yout] = ode45(@(t,x) model9(t,x,am,bm,gs,as,bs,ap,bp,k1,k2),[0 TLim],x0); 

figure(4) 
set(gca,'FontSize',35);
plot(Tout,Yout(:,1),tDat,mDat,'LineWidth',5);
xlabel('Time (hours)'); ylabel('Count'); 
legend('m Sim','m Real');

figure(5) 
set(gca,'FontSize',35);
plot(Tout,Yout(:,2),tDat,sDat,'LineWidth',5);
xlabel('Time (hours)'); ylabel('Count'); 
legend('s Sim','s Real');

figure(6) 
set(gca,'FontSize',35);
plot(Tout,Yout(:,3),tDat,pDat,'LineWidth',5);
xlabel('Time (hours)'); ylabel('RFU'); 
legend('p Sim','p Real');


saveas(4,strcat('.\',SavDir,'\model6mRNA.jpg'));
saveas(5,strcat('.\',SavDir,'\model6miRNA.jpg'));
saveas(6,strcat('.\',SavDir,'\model6protein.jpg'));

save(strcat('.\',SavDir,'\ModelFitIFFL9'),'am','bm','gs','k1','k2','ap','bp','as','bs');