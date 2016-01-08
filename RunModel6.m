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

mRange = 'D74:D79';
sRange = 'D80:D85'; 
tRange = 'C80:C85';
pRange = 'F86:F91';
datFile = 'AlexIFFLdata.xls';
SavDir = 'Model6d123'; 
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
yMin = [100 .01 0 0 0 0 .001]; 
yMax = [3000 .5 .001 .9999 100 5000 1];

y0 = [248 .19 3.2256e-4 .01 .1 45 .03]; 

y = fmincon(@(y) model6mp(y,sDat,mDat,pDat,tDat,as,bs),y0,[],[],[],[],yMin,yMax,[],options); 

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

[Tout,Yout] = ode45(@(t,x) model6(t,x,am,bm,gs,as,bs,ap,bp,k1,k2),[0 TLim],x0); 

figure(4) 
plot(Tout,Yout(:,1),tDat,mDat,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('m Sim','m Real');

figure(5) 
plot(Tout,Yout(:,2),tDat,sDat,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('s Sim','s Real');

figure(6) 
plot(Tout,Yout(:,3),tDat,pDat,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('p Sim','p Real');


saveas(4,strcat('.\',SavDir,'\model6mRNA.jpg'));
saveas(5,strcat('.\',SavDir,'\model6miRNA.jpg'));
saveas(6,strcat('.\',SavDir,'\model6protein.jpg'));

save(strcat('.\',SavDir,'\ModelFitIFFL6'),'am','bm','gs','k1','k2','ap','bp','as','bs');