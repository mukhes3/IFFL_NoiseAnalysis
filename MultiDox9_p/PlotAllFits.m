%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Plot all fits in the same figure 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear; 
close all; 

%%                       Run Code 

%666 

load('666model4params'); 

[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 72],x0); 
p666=y(:,3); 
t666 = t;
p666a = pDat; 
t666a = tDat; 

%800 

load('800model4params'); 

[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 72],x0); 
p800=y(:,3); 
t800 = t;
p800a = pDat; 
t800a = tDat; 

%1000 

load('1000model4params'); 

[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 72],x0); 
p1000=y(:,3); 
t1000 = t;
p1000a = pDat; 
t1000a = tDat; 

%1500 

load('1500model4params'); 

[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 72],x0); 
p1500=y(:,3); 
t1500 = t;
p1500a = pDat; 
t1500a = tDat; 

%2000 

load('2000model4params'); 

[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 72],x0); 
p2000=y(:,3); 
t2000 = t;
p2000a = pDat; 
t2000a = tDat; 

%2500 

load('2500model4params'); 

[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 72],x0); 
p2500=y(:,3); 
t2500 = t;
p2500a = pDat; 
t2500a = tDat; 

%5000 

load('5000model4params'); 

[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 72],x0); 
p5000=y(:,3); 
t5000 = t;
p5000a = pDat; 
t5000a = tDat; 

%15000 

load('15000model4params'); 

[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 72],x0); 
p15000=y(:,3); 
t15000 = t;
p15000a = pDat; 
t15000a = tDat; 

%1000000 

load('1000000model4params'); 

[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 72],x0); 
p1000000=y(:,3); 
t1000000 = t;
p1000000a = pDat; 
t1000000a = tDat; 


%%              Plotting 


figure(1)
set(gca,'FontSize',35);
hold on; 
plot(t666,p666,t800,p800,t1000,p1000,t1500,p1500,t2000,p2000,t2500,p2500,t5000,p5000,t15000,p15000,t1000000,p1000000,t666a,p666a,'k-.',t800a,p800a,'k-.',t1000a,p1000a,'k-.',t1500a,p1500a,'k-.',t2000a,p2000a,'k-.',t2500a,p2500a,'k-.',t5000a,p5000a,'k-.',t15000a,p15000a,'k-.',t1000000a,p1000000a,'k-.','LineWidth',3); 
xlabel('Time (hours)'); ylabel('Fluorescence (RFU)'); 
% legend('Dox-666','Dox-800','Dox-1000','Dox-1500','Dox-2000','Dox-2500','Dox-5000','Dox-15000','Dox-1000000','Dox-666 actual','Dox-800 actual','Dox-1000 actual','Dox-1500 actual','Dox-2000 actual','Dox-2500 actual','Dox-5000 actual','Dox-15000 actual','Dox-1000000 actual','FontSize',20,'Location','BestOutside'); 
saveas(1,'MultiDoxPlotActualIncl.jpg');
