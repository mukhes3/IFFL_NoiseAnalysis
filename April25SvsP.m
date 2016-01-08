%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Synthetic S vs P plot 
%                                   Sumit Mukherjee 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear; 
close all; 

%%                          Initialization 

LoadFile = 'IFFL_Alex_pv_FISH_mod\model4params'; 
load(LoadFile); 
SavFile = '.\IFFL_Alex_pv_FISH_mod\SvsP.jpg'; 

S = logspace(-3,5,50); 

P = (ap*am/bp)./(bm + gs*S);

I = find(P<= .5*(ap*am/(bp*bm))); 
I50 = S(I(1)); 

figure(1) 
semilogx(S,P,'LineWidth',3)
title(strcat('I_{50} =',num2str(I50,'%4.4f\n')));
set(gca,'FontSize',15);
saveas(1,SavFile); 