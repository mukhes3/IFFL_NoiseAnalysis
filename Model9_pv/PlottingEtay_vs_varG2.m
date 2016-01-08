%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program generates a figure of how eta_y varies with variance 
%   in Var(g_2^ss) for different values of gamma_s 
%
%                                             Sumit Mukherjee
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear; 
close all; 

%%                           Initialization 

load('ModelFitIFFL9'); 
K = k1*k2; 
d = as/bs; 

Eg2 = .9; 


%%                              Plotting

Varg2 = [.01:.01:.8]; 
GammaS = [1e-3]; 
nY = zeros(length(GammaS),length(Varg2)); 

for i = 1:length(GammaS)
    gs = GammaS(i); 
    Nu = (((gs*d - K*d*bm)/(1 + K*d*Eg2)^2)^2)*Varg2; 
    Den = ((bm + gs*d*Eg2)/(1 + K*d*Eg2)  - (K*d^2*(gs-K*bm)/(1 + K*d*Eg2)^3)*Varg2).^2; 
    nY(i,:) = Nu./Den; 
end

plot(Varg2,nY)


