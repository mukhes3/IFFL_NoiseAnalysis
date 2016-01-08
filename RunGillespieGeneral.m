%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%             Generalized program for Gillespie simulation
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
%%             User Specified Compoments 

% Input Stoichiometric Matrix 

Si = [ 1 1 0 ; ...
       0 1 1 ; ... 
       1 0 1 ]; 
   
% Output Stoichiometric Matrix 

So = [ 2 0 0 ; ...
       0 2 0 ; ... 
       0 0 2 ];
   
% Rate Constant Matrix 

k = 10*[ 1 1 1 ]'; 

% Initial Number of molecules of species 

x0 = [ 5 10 15 ]';

% Run time 

Tmax = 1e3 ; 

 
%%            Run Simulation 

[Xrec, trec] = GillespieGeneral(Si, So, k, x0, Tmax) 

%%           Generate Plots 

figure(1) 
plot(trec,Xrec,'LineWidth',3);
xlabel('time');
ylabel('Number of Molecules');

saveas(1,'GillespieGeneral.jpg');










