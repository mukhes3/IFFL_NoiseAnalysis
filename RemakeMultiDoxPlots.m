%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          This program plots the multiple Dox level plots 
%
%                                               Sumit Mukherjee
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear; 
close all; 

%%                       Initialization 

datFile = 'OthersMultiDox.xls';
r400 = 'D1:D9'; 
r666 = 'D10:D18'; 
r800 = 'D19:D27'; 
r1000 = 'D28:D36'; 
r1500 = 'D37:D45'; 
r2000 = 'D46:D54';
r2500 = 'D55:D63'; 
r5000 = 'D64:D72'; 
r15000 = 'D73:D81'; 
r1000000 = 'D82:D90'; 


t400 = 'C1:C9'; 
t666 = 'C10:C18'; 
t800 = 'C19:C27'; 
t1000 = 'C28:C36'; 
t1500 = 'C37:C45'; 
t2000 = 'C46:C54';
t2500 = 'C55:C63'; 
t5000 = 'C64:C72'; 
t15000 = 'C73:C81'; 
t1000000 = 'C82:C90'; 

matFile = '.\IFFL_Alex_d23\model4params'; 

OutFolderName = 'MultiDox_d23'; 

%%               Re-make plots 


%666
OutFileAdd = '666';
tRange = t666;
pRange = r666;
run('RunMultiDox.m'); 

%800
OutFileAdd = '800';
tRange = t800;
pRange = r800;
run('RunMultiDox.m'); 

%1000
OutFileAdd = '1000';
tRange = t1000;
pRange = r1000;
run('RunMultiDox.m'); 

%1500
OutFileAdd = '1500';
tRange = t1500;
pRange = r1500;
run('RunMultiDox.m'); 

%2000
OutFileAdd = '2000';
tRange = t2000;
pRange = r2000;
run('RunMultiDox.m'); 

%2500
OutFileAdd = '2500';
tRange = t2500;
pRange = r2500;
run('RunMultiDox.m'); 

%5000
OutFileAdd = '5000';
tRange = t5000;
pRange = r5000;
run('RunMultiDox.m'); 

%15000
OutFileAdd = '15000';
tRange = t15000;
pRange = r15000;
run('RunMultiDox.m'); 

%1000000
OutFileAdd = '1000000';
tRange = t1000000;
pRange = r1000000;
run('RunMultiDox.m'); 