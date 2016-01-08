%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Surface Plot Protein levels 
%                                           Sumit Mukherjee
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear; 
close all; 

%%                     Run Code 

load('Pparams'); 
Z0 = pVal; 

load('PVparams'); 
Z4 = pVal; 

load('d3params'); 
Z3 = pVal; 

load('d23params'); 
Z2 = pVal;

Z = [Z0;Z2;Z3;Z4]; 

X = [0 2 3 4]; 
Y = Dox; 

surf(X,Y,Z')

