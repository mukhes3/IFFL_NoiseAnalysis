%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Paramter Estimation of mRNA from control case
%                                                   Sumit Mukherjee
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clc; 
 clear;
 close all; 
 
 %%                      Initialization

mRange = 'D20:D25';
sRange = 'D26:D31'; 
tRange = 'C20:C25';
datFile = 'AlexIFFLdata.xls';


mDat = xlsread(datFile,mRange); 
sDat = xlsread(datFile,sRange);
tDat = xlsread(datFile,tRange);

%%                  Parameter estimation for mRNA 

yMin = [0 .05]; 
yMax = [500 .5];

y0 = [40 .10]; 

options=optimset('Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolCon',1e-9,'TolFun',1e-9,'TolX',1e-9);

y = fmincon(@(y) model2sdat(y,mDat,tDat),y0,[],[],[],[],yMin,yMax,[],options); 

am = y(1); 
bm = y(2); 


% am = 193.5477; 
% bm = .1640;
% gs = .0003; 

%%               Plotting 

m = mDat(1);
mRec = m; 

for i=1:5

T = tDat(i+1)-tDat(i);    

m = am*T + (1-bm*T)*m;  

mRec = [mRec;m];

end


figure(1) 
plot(tDat,mDat,tDat,mRec,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('m Real','m Sim');

saveas(1,'model3mRNA.jpg');


