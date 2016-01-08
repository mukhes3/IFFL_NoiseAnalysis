%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Noise Rejection Analysis 
%                                  Sumit Mukherjee
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear; 
close all; 

%%                          initialization 

LoadFile = 'IFFL_Alex_pv_FISH_mod\model4params'; 
simFile = 'McherryRun'; 
savFolder = 'BufferObserve'; 
load(LoadFile); 

% Realistically chosen initial conditions
m0 = 100; 
s0 = 80; 
p0 = 1800; 

%%         Calculating Atenuation at different frequencies for IFFL


T = fliplr(logspace(-5,1,100)); 
F = 2*pi./T; 

N = length(T); 
Arec = zeros(N,1); 
Aprec = zeros(N,1);

for i = 1:N 
    f = F(i); 
    Tmax = max(1000*T(i),300); 
    Tsamp = T(i)/25;
    temp = sim('McherryRun'); 
    l = length(mRec.time); 
    In1 = ceil(.9*l); 
    A = .5*abs(max(mRec.signals.values(In1:l)) - min(mRec.signals.values(In1:l))); 
    Ap = .5*abs(max(pRec.signals.values(In1:l)) - min(pRec.signals.values(In1:l)));
    Aprec(i) = Ap;
    Arec(i) = A; 
end

%%       Calculating Atenuation at different frequencies for Open Loop

gs = 0; 

Arec2 = zeros(N,1); 
Aprec2 = zeros(N,1);
for i = 1:N 
    f = F(i); 
    Tmax = max(1000*T(i),300); 
    Tsamp = T(i)/25;
    temp = sim('McherryRun'); 
    l = length(mRec.time); 
    In1 = ceil(.9*l); 
    A = .5*abs(max(mRec.signals.values(In1:l)) - min(mRec.signals.values(In1:l))); 
    Ap = .5*abs(max(pRec.signals.values(In1:l)) - min(pRec.signals.values(In1:l)));
    Aprec2(i) = Ap;
    Arec2(i) = A; 
end

%%      Comparing atenuations 

figure(1) 
semilogx(F,20*log10(Arec),F,20*log10(Arec2),'LineWidth',3); 
xlabel('Frequency'); ylabel('Amplitude of oscillation'); 
title('Open Loop vs IFFL attenuation'); 
legend('IFFL','Open Loop','Location','BestOutside'); 

saveas(1,strcat('.\',savFolder,'\OLvsIFFLattn.jpg')); 

figure(2) 
semilogx(F,20*log10(Aprec),F,20*log10(Aprec2),'LineWidth',3); 
xlabel('Frequency'); ylabel('Amplitude of oscillation'); 
title('Open Loop vs IFFL attenuation Protein'); 
legend('IFFL','Open Loop','Location','BestOutside'); 

saveas(2,strcat('.\',savFolder,'\OLvsIFFLattnProtein.jpg')); 
