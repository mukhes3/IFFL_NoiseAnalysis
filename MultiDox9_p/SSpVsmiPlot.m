%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Plot Steady State Protein VS miRNA levels
%                                                   Sumit Mukherjee
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear; 
close all; 

%%                      Reading Values 

load('666model4params'); 

p666a = pDat(end); 
s666a = as/bs; 

%800 

load('800model4params'); 

p800a = pDat(end); 
s800a = as/bs; 

%1000 

load('1000model4params'); 

p1000a = pDat(end); 
s1000a = as/bs;  

%1500 

load('1500model4params'); 

p1500a = pDat(end); 
s1500a = as/bs; 
%2000 

load('2000model4params'); 

p2000a = pDat(end); 
s2000a = as/bs; 

%2500 

load('2500model4params'); 

p2500a = pDat(end); 
s2500a = as/bs; 

%5000 

load('5000model4params'); 

p5000a = pDat(end); 
s5000a = as/bs; 

%15000 

load('15000model4params'); 

p15000a = pDat(end); 
s15000a = as/bs; 

%1000000 

load('1000000model4params'); 

p1000000a = pDat(end); 
s1000000a = as/bs; 

%%                               Plotting 


p = [p666a p800a p1000a p1500a p2000a p2500a p5000a p15000a p1000000a]; 

s = [s666a s800a s1000a s1500a s2000a s2500a s5000a s15000a s1000000a];

plot(s,p)