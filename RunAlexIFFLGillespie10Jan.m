%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%             Generalized program for Gillespie simulation
%                                                       Sumit Mukherjee 
% 
%         Things to change when you run the code 
%  1) RunFileName 
%  2) OutFolderName - Name of the folder where the image will be stored 
%  3) OutFileAdd - prefix of the 3 images (the suffix is safe for all i.e.
%  mRNA, miRNA and protein.jpg) 
%  4) k - This contains the names of the variables. The need to update this
%  variable might be eliminated in future versions 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all; 
%%             User Specified Compoments 

ReadFileName = '.\IFFL_Alex_pv_FISH\model4paramsUpdated2'; 
OutFolderName = 'Gillespie_Alex_pv_FISH'; 
OutFileAdd = 'Updated2';


% Input Stoichiometric Matrix 

Si = [ 0 0 0 ; ...
       0 0 0 ; ... 
       1 0 0 ; ...
       1 1 0 ; ...
       0 1 0 ; ...
       1 0 0 ; ...
       0 0 1 ...
       ]; 
   
% Output Stoichiometric Matrix 

So = [ 1 0 0 ; ...
       0 1 0 ; ... 
       0 0 0 ; ...
       0 1 0 ; ...
       0 0 0 ; ...
       1 0 1 ; ...
       0 0 0 ...
       ]; 

   load(ReadFileName); 

% Rate Constant Matrix (this has to be changed according to names of the
% variables used) 

k = [ amo2 as bmo2 gso2 bs apo2 bpo2 ]'; 

% Initial Number of molecules of species 

x0 = [ 5 10 15 ]';

% Run time 

Tmax = 100 ; 



%%            Run Simulation 

[Xrec, trec] = GillespieGeneral(Si, So, k, x0, Tmax); 

%%           Generate Plots 

figure(1) 
plot(trec,Xrec(1,:),'LineWidth',3);
xlabel('time');
ylabel('mRNA molecules');

figure(2) 
plot(trec,Xrec(2,:),'LineWidth',3);
xlabel('time');
ylabel('miRNA molecules');

figure(3) 
plot(trec,Xrec(3,:),'LineWidth',3);
xlabel('time');
ylabel('protein molecules');

mkdir(OutFolderName);
saveas(1,strcat(OutFolderName,'\',OutFileAdd,'mRNA.jpg'));
saveas(2,strcat(OutFolderName,'\',OutFileAdd,'miRNA.jpg'));
saveas(3,strcat(OutFolderName,'\',OutFileAdd,'protein.jpg'));










