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

%%                     Parameter Estimation  for miRNA

xMin = [0 .02]; 
xMax = [500 .04]; 
% xMax = 1e3*ones(4,1); 
% phi = [am bm gs as bs];

x0 = [ 45 .03]; 

options=optimset('Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolCon',1e-9,'TolFun',1e-9,'TolX',1e-9);

x = fmincon(@(x) model2sdat(x,sDat,tDat),x0,[],[],[],[],xMin,xMax,[],options); 

 
as = x(1); 
bs = x(2); 

%%                 gs Parameter estimation for mRNA 

yMin = [0]; 
yMax = [1];

y0 = [.0001]; 

am = 142.3237; 
bm = .0522;

y = fmincon(@(y) model3gs(y,mDat,sDat,tDat,as,bs,am,bm),y0,[],[],[],[],yMin,yMax,[],options); 

gs = y(1); 

% am = 193.5477; 
% bm = .1640;
% gs = .0003; 
% 
% gs = .00; 
%%               Plotting 

s = sDat(1); 
m = mDat(1);
sRec = s; 
mRec = m; 

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

saveas(1,'model3mRNA.jpg');



figure(2)
plot(tDat,sDat,tDat,sRec,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('s Real','s Sim');

saveas(2,'model3miRNA.jpg');