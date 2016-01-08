a = random('poiss',1*ones(100000,1)); 
b = 1 + 100*random('poiss',100*ones(100000,1)); 

c=[]; 

% for i = 1:100
%     for j = 1:100
%         if b(j)>0
%             c = [c,a(i)/b(j)];
%          
%         end
%     end
% end
c = a./b; 


act_cv = std(a)/mean(a); 
cv_c = std(c)/mean(c); 