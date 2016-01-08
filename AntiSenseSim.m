%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Model Anti Sense Addition 
%                                              Sumit Mukherjee 
%
% t1 = time at which the antisense miRNA's are transfected (they are
% assumed to cause a alpha times reduction in the miRNA numbers). We also
% assume that the concentration of the antisense miRNA's is high enough
% that their reduction over time is negligible. 
%
% t2 = time at which mCherry gene is induced. Once again 0th order dynamics
% are assumed here. i.e. the 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear; 
close all; 

%%                               Initialization 

ParamFile = 'IFFL_Alex_pv_March5modFISh\model4params'; 
load(ParamFile); 

t1 = 0; 
t2 = 60; 

Tmax = 250;

m0 = 63; 
s0 = 105; 
p0 = 2700; 

x0 = [m0; s0; p0]; 
alpha = 0.9; 
%%                Running ODE 

asd = (1-alpha)*as; 
if t1 <= t2 
    [Tout,Yout] = ode45(@(t,x) model2sdatcont(t,x,am,bm,gs,asd,bs,ap,bp),[0 Tmax-t1],x0); 
else 
    [Tout1,Yout1] = ode45(@(t,x) model2sdatcont(t,x,am,bm,gs,as,bs,ap,bp),[0 t1],x0);
    x0 = Yout1(end,:)';
    [Tout2,Yout2] = ode45(@(t,x) model2sdatcont(t,x,am,bm,gs,asd,bs,ap,bp),[0 Tmax-t1],x0);
    Tout = [Tout1; (Tout2 + t1)]; 
    Yout = [Yout1;Yout2]; 
end

%%                Plotting 

if t1<=t2 
    
figure(1) 
plot(Tout,Yout(:,1),'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
title('mRNA vs Time'); 

figure(2) 
plot(Tout,Yout(:,2),'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
title('miRNA vs Time'); 

figure(3) 
plot(Tout,Yout(:,3),'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
title('Protein vs Time'); 

else 
    
figure(1) 
plot(Tout1,Yout1(:,1),(Tout2+t1),Yout2(:,1),'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('Pre Addition','Post Addition','Location','BestOutside'); 
title('mRNA vs Time'); 

figure(2) 
plot(Tout1,Yout1(:,2),(Tout2+t1),Yout2(:,2),'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.');
legend('Pre Addition','Post Addition','Location','BestOutside');
title('miRNA vs Time'); 

figure(3) 
plot(Tout1,Yout1(:,3),(Tout2+t1),Yout2(:,3),'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('Pre Addition','Post Addition','Location','BestOutside');
title('Protein vs Time'); 
end