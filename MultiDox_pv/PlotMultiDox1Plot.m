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
p666=y(end,3); 
figure(1) 
hold on; 
plot(t,y(:,3),'LineWidth',3,'Color','b'); 
xlabel('Time (hours)'); ylabel('Protein Conc.'); 


load('Dox800model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
plot(t,y(:,3),'LineWidth',3,'Color','g'); 
p800=y(end,3); 

load('Dox1000model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
plot(t,y(:,3),'LineWidth',3,'Color','r'); 
p1000=y(end,3); 


load('Dox1500model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
plot(t,y(:,3),'LineWidth',3,'Color','c'); 
p1500=y(end,3); 

load('Dox15000model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
plot(t,y(:,3),'LineWidth',3,'Color','m'); 
p15000=y(end,3); 

load('Dox1000000model4params'); 
x0 = [63; 105; 2700]; 
[t,y] = ode45(@(t,x) mpdat4(t,x,am,bm,gs,ap,bp,as,bs), [0 150],x0); 
plot(t,y(:,3),'LineWidth',3,'Color','y'); 
p1000000=y(end,3); 

legend('Dox-666','Dox-800','Dox-1000','Dox-1500','Dox-15000','Dox-1000000','Location','BestOutside'); 
saveas(1,'MultiDoxPlot.jpg');

%%          Plot Steady State Protein Level Vs Dox level 

figure(2)
Dox = [666 800 1000 1500 15000 1000000]; 
pVal = [p666 p800 p1000 p1500 p15000 p1000000]; 
semilogx(Dox,pVal,'LineWidth',3); 
xlabel('Dox level in ng/ul'); ylabel('Fluorescence (RFU)'); 
saveas(2,'DoxVsProteinIFFL.jpg'); 
