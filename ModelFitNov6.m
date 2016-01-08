%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%             Model Fitting IFFL 
%                                  Sumit Mukherjee 
%
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear; 
close all; 

%%                      Initialization

mRange = 'D2:D7';
sRange = 'D8:D13'; 
tRange = 'C2:C7';
datFile = 'AlexIFFLdata.xls';


mDat = xlsread(datFile,mRange); 
sDat = xlsread(datFile,sRange);
tDat = xlsread(datFile,tRange);

%%                     Parameter Estimation  

xMin = zeros(5,1); 
xMax = [500 10 10 500 10]; 
% xMax = 1e3*ones(4,1); 
% phi = [am bm gs as bs];

x0 = [ 50 .1 .00075 .4*50 .03 ]; 

options=optimset('Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolCon',1e-9,'TolFun',1e-9,'TolX',1e-9);

x = fmincon(@(x) model1(x,mDat,sDat,tDat),x0,[],[],[],[],xMin,xMax,[],options); 

am = x(1); 
bm = x(2); 
gs = x(3); 
as = x(4); 
bs = x(5); 

%%               Plotting 

m = mDat(1); 
s = sDat(1); 

mRec = m;
sRec = s; 

for i=1:5

T = tDat(i+1)-tDat(i);    

m = am*T + (1-bm*T)*m - gs*m*s*T; 
s = as*T + (1-bs*T)*s ;

mRec = [mRec;m];
sRec = [sRec;s];

end


figure(1) 
plot(tDat,mDat,tDat,mRec,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('m Real','m Sim');

figure(2) 
plot(tDat,sDat,tDat,sRec,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('s Real','s Sim');

saveas(1,'model1mRNA.jpg');
saveas(2,'model1miRNA.jpg');

save('ModelNov6Dat','am','bm','gs','as','bs'); 
