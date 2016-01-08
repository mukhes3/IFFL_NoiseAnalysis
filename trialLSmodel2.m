clc;
clear; 

mRange = 'D2:D7';
sRange = 'D8:D13';
tRange = 'C2:C7';
datFile = 'AlexIFFLdata.xls';
mDat = xlsread(datFile,mRange);
sDat = xlsread(datFile,sRange);
tDat = xlsread(datFile,tRange);
D = mDat(3:end)-mDat(2:5)
T = tDat(3:6)-tDat(2:5)
m = mDat(2:5)
s = sDat(2:5)
A = [T, -m.*T, -m.*s.*T]
inv(A'*A)*A'*D