%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Theoretical Validation (from SHS) of simulation results
%                                                  Sumit Mukherjee
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear; 
close all; 

%%               Initialization 

LoadFileName ='.\IFFL_Alex_pv_FISH\model4params'; 
load(LoadFileName)


%%               Calculation 

% Ss = as/bs
% 
% Ms = am/(bm + gs*Ss)

% open 

lp = .8; 
lm = .18;
Gmax = 1;
K = lp+lm; 
T1 = inv(lp+lm); 
T2 = inv(bm); 
Pon = lp*T1; 

%FromGeorg
MsGeorg = (lp/(lp + lm))*am/bm
etaM2Georg = (1/MsGeorg) + (lm/lp)*(bm/(bm + K))
etaMGeorg = sqrt(etaM2Georg)

%From Paulsson 
GsPaul = lp*Gmax*T1
MsPaul = T2*am*GsPaul 

etaGPaul = sqrt((1-Pon)/GsPaul)
etaMPaul = sqrt((1/MsPaul) + ((1-Pon)/GsPaul)*(T1/(T1+T2)))

%%              Closed loop 

Gc = lp*Gmax/K
Sc = as*Gc/bs 
Tm = bm + gs*Sc 
q = gs*Sc/Tm 
Mc = am*Gc/Tm 

EtaM2 = inv(Mc*(1-q)) + (lm/lp)*(bm/(K+bm))*( 1 - q - (bm*K/((K+bs)*(bm+bs)))*q*((2-q)/(1-q))) + inv(Sc*(bm+bs)*(1-q))*bm*q^2

