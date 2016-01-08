%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Program to plot different protein plots for different Dox levels 
%                                                 Sumit Mukherjee 
% Note :- This is a little hardcoded 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear; 
close all; 

load('Dox666model4params'); 

x0 = [63; 105; 2700]; 

[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
p666=y(:,3); 
m666=y(:,1);
t666 = t;


load('Dox800model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
p800=y(:,3); 
m800=y(:,1);
t800=t;

load('Dox1000model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
p1000=y(:,3); 
m1000=y(:,1);
t1000=t; 

load('Dox1500model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
p1500=y(:,3); 
m1500=y(:,1);
t1500 = t; 

load('Dox2000model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
p2000=y(:,3); 
m2000=y(:,1);
t2000=t; 

load('Dox2500model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
p2500=y(:,3); 
m2500=y(:,1);
t2500 = t; 


load('Dox5000model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
p5000=y(:,3); 
m5000=y(:,1);
t5000 = t; 

load('Dox15000model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
p15000=y(:,3); 
m15000=y(:,1);
t15000=t; 

load('Dox1000000model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
p1000000=y(:,3); 
m1000000=y(:,1);
t1000000=t; 

figure(1)
set(gca,'FontSize',35);
plot(t666,p666,t800,p800,t1000,p1000,t1500,p1500,t2000,p2000,t2500,p2500,t5000,p5000,t15000,p15000,t1000000,p1000000,'LineWidth',3); 
xlabel('Time (hours)'); ylabel('Fluorescence (RFU)'); 
legend('Dox-666','Dox-800','Dox-1000','Dox-1500','Dox-2000','Dox-2500','Dox-5000','Dox-15000','Dox-1000000','Location','BestOutside'); 
saveas(1,'MultiDoxPlot.jpg');

%%          Plot Steady State Protein Level Vs Dox level 

figure(2)
Dox = [666 800 1000 1500 2000 2500 5000 15000 1000000]; 
pVal = [p666(end) p800(end) p1000(end) p1500(end) p2000(end) p2500(end) p5000(end) p15000(end) p1000000(end)]; 
semilogx(Dox,pVal,'LineWidth',3); 
xlabel('Dox level in ng/ul'); ylabel('Fluorescence (RFU)'); 
saveas(2,'DoxVsProteinIFFL.jpg'); 
save('e:\UWPhDWork\ModelFittingIFFL\PVparams','Dox','pVal'); 


figure(3)
plot(t666,m666,t800,m800,t1000,m1000,t1500,m1500,t2000,m2000,t2500,m2500,t5000,m5000,t15000,m15000,t1000000,m1000000,'LineWidth',3); 
xlabel('Time (hours)'); ylabel('mRNA no.'); 
legend('Dox-666','Dox-800','Dox-1000','Dox-1500','Dox-2000','Dox-2500','Dox-5000','Dox-15000','Dox-1000000','Location','BestOutside'); 
saveas(3,'MultiDoxPlotmRNA.jpg');

figure(4) 