load('.\IFFL_Alex_d3\model4params'); 


TLim = 100; 

x0 = [60; 70; 1000]; 

[Tout1,Yout1] = ode45(@(t,x) model2sdatcont(t,x,am,bm,gs,as,bs,ap,bp),[0 TLim],x0); 
[Tout2,Yout2] = ode45(@(t,x) model2sdatcont(t,x,am,bm,0,as,bs,ap,bp),[0 TLim],x0);


figure(1) 
plot(Tout1,Yout1(:,3),Tout2,Yout2(:,3),'LineWidth',3);
xlabel('Time (hours)'); ylabel('Conc.'); 
legend('Before RBP knockdown','After RBP knockdown','Location','Best');


saveas(1,'APAregDemo.jpg');

