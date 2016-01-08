%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             How aplha m varies with induction levels 
%                                                Sumit Mukherjee
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear; 
close all;

%%

DoxLevels = [666 800 1000 1500 2000 2500 5000 15000 1000000]; 
amRec = zeros(1,length(DoxLevels)); 

for i = 1:length(DoxLevels)
    load(strcat(int2str(DoxLevels(i)),'model4params')); 
    amRec(i) = am ; 
end

figure(1) 

semilogx(DoxLevels,amRec,'LineWidth',3);
xlabel('Dox level (in ng)','FontSize',25); 
ylabel('\alpha_m','FontSize',25);
title('\alpha_m vs Dox level plot for PV','FontSize',25); 
set(gca,'FontSize',15); 
