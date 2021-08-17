function showAcc()

[u2,~,~,~,~,~]=simpleAP_Matrix(120,800000,1,25);
%[rows2,cols2]=size(u2)
nT = [20000:40000:780000];
for j =1:length(nT)
%[u1,~,~,~,~,~]=simpleAP(120,nT(j),1,25);

[u1,~,~]=opSpltC(120,nT(j),0,1,0,25);
%[~,cols1]=size(u1)
%mid1 = (cols1-1)/2;
%mid2 = (cols2-1)/2;
size(u1(:,end))
size(u2(:,end))
acc = abs(u1(1:121,end)-u2(:,end));
maxAcc(j)=max(acc);
end

loglog(1./nT,maxAcc,'Linewidth',2)
xlabel('dt')
ylabel('acc')
title('Accuracy','fontsize',24);
set(gca,'FontSize',18)
% plot(x1,maxAcc,'r','Linewidth',2);
% xlabel('x pos','fontsize',24);
%  ylabel('Diff','fontsize',24);
%  title('Max Difference between both MOL schemes','fontsize',24);
%for j=1:cols
%    plot(x1,acc(:,j));
%    drawnow
%end
end

