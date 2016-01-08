OpReads = xlsread('FISHAlex.xls','A258:A294'); 
IfflReads1 = xlsread('FISHAlex.xls','A629:A674');
IfflReads2 = xlsread('FISHAlex.xls','A717:A765');
IfflReads = [IfflReads1;IfflReads1];
mean(OpReads)
mean(IfflReads)
mean(IfflReads1)

figure(1)
set(gca,'FontSize',20);
hist(OpReads,20); 
xlabel('mRNA reads'); 
ylabel('No. of instances'); 
saveas(1,'histFishOp.jpg')

figure(2)
set(gca,'FontSize',20);
hist(IfflReads,20); 
xlabel('mRNA reads'); 
ylabel('No. of instances'); 
saveas(2,'histFishIFFL.jpg')




