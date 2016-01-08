%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%             Model Fitting IFFL 
%                                  Sumit Mukherjee 
%
%  In this code the miRNA and proten data's are fit separately and then the 
% mRNA data is fit by choosing several different initial conditions and
% then picking the minimum of the lot. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear; 
close all; 

%%                      Initialization

mRange = 'D2:D7';
sRange = 'D8:D13'; 
tRange = 'C2:C7';
pRange = 'D14:D19';
datFile = 'AlexIFFLdata.xls';


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

options=optimset('Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolCon',1e-9,'TolFun',1e-9,'TolX',1e-9);

x = fmincon(@(x) model2sdat(x,sDat,tDat),x0,[],[],[],[],xMin,xMax,[],options); 

 
as = x(1); 
bs = x(2); 

%%                   Parameter Estimation for mCherry protein

xMin = [0 .001]; 
xMax = [5000 1]; 
% xMax = 1e3*ones(4,1); 
% phi = [am bm gs as bs];

x0 = [ 45 .03]; 

options=optimset('Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolCon',1e-9,'TolFun',1e-9,'TolX',1e-9);

x = fmincon(@(x) model2pdat(x,mDat,pDat,tDat),x0,[],[],[],[],xMin,xMax,[],options); 

 
ap = x(1); 
bp = x(2); 



%%                  Parameter estimation for mRNA 

yMin = [0 .01 0]; 
yMax = [3000 .5 .001];

yRec=zeros(100,3); 
fvalRec = zeros(100,1); 

for i = 1:100 
    
y0 = yMin + (yMax-yMin)*rand(3,1); 

[yRec(i,:),fvalRec(i)] = fmincon(@(y) model2mdat(y,mDat,sDat,tDat,as,bs),y0,[],[],[],[],yMin,yMax,[],options); 

end

[fval,index] = min(fvalRec); 

am = yRec(index,1); 
bm = yRec(index,2); 
gs = yRec(index,3); 

% am = 193.5477; 
% bm = .1640;
% gs = .0003; 

%%               Plotting 

s = sDat(1); 
m = mDat(1);
p = pDat(1);
sRec = s; 
mRec = m; 
pRec = p; 

for i=1:5

T = tDat(i+1)-tDat(i);    

m = am*T + (1-bm*T)*m - gs*m*s*T;  
s = as*T + (1-bs*T)*s ;
p = ap*T*mDat(i) + (1-bp*T)*p ;

mRec = [mRec;m];
sRec = [sRec;s];
pRec = [pRec;p];

end


figure(1) 
plot(tDat,mDat,tDat,mRec,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('m Real','m Sim');

saveas(1,'model2mRNA.jpg');



figure(2)
plot(tDat,sDat,tDat,sRec,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('s Real','s Sim');

saveas(2,'model2miRNA.jpg');

figure(3) 
plot(tDat,pDat,tDat,pRec,'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('p Real','p Sim');

saveas(3,'model2p.jpg');