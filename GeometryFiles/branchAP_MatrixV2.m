function [u,mm,nn,hh,x,t]=branchAP_MatrixV2(nX,nT,X,T)
% nX = number of spatial steps along a rod
% nT = number of time steps
% X = Length of rod
% T = end time

%Test case: [u,mm,nn,hh,x,t]=simpleAP(100,200000,1,20);
%Test case: [u,mm,nn,hh,x,t]=simpleAP(200,200000,2,20);
%Test case: [u,mm,nn,hh,x,t]=simpleAP(200,200000,2,30);
if nargin == 0
    [nT, X, T]=deal(40000,1,100);
    %loc = 59;
    src_loc=1;
end

[~,Inter_node,End_node]=get_IntrsctNode();

x1=linspace(0,1,End_node);
x2=x1(Inter_node+1:end);
size(x1);
size(x2);

x=[x1,x2];

 k=T/nT; % time step size
% h=X/nX; % equidistant spatial discretization width
% x=linspace(0,X,nX+1);    % space vector
 t=linspace(0,T,nT+1);    % time vector
 
 h = X/(End_node-1);
 nX=length(x);

%----------------------------------------%
%Parameters for AP
a=.001;      %radius
R=10;       %resistance
gk=36;      %potassium ion conductance
ek=-12;     %potassium reversal potential
gna=153;  %120  %sodium ion conductance
ena=112;    %115; %sodium reversal potential
ni=0.5; 
mi=0.4;
hi=0.2;
gl=0.3;     %leak conductance
el=0.6;    %leak potential

%compute coefficients
b=a/(2*R);
c=0.09;%c=0.9e-2;
%----------------------------------------%
fprintf('Number of spatial steps = %i\n',nX);
fprintf('Spatial dx = %d\n',h);
fprintf('Number of time steps = %i\n',nT);
fprintf('Temporal dt = %d\n',k);
% 
%Initialize Solutions
u=zeros(nX,nT+1);
nn=zeros(nX,nT+1);
mm=zeros(nX,nT+1);
hh=zeros(nX,nT+1);
%set boundary conditions at one end of line;
vstart=-55;
%u(1,1:ceil((nT+1)/2))=-vstart;
u(src_loc,1:end)=-vstart*0;
% 
nn(1,1)=ni;
mm(1,1)=mi;
hh(1,1)=hi;
% %Method of Lines Computation
% %The i-th column is a time slice
% 
% %define stencil matrix
sten = [1 -2 1];
sysMat = spdiags(ones(End_node,1)*sten,-1:1,End_node,End_node)*b/h^2; %sparse diagonal
size(x1);
size(x2);
size(sysMat);
rem = nX-End_node;
sysMat2 = spdiags(ones(rem,1)*sten,-1:1,rem,rem)*b/h^2; %sparse diagonal
size(sysMat2);
Mat=blkdiag(sysMat,sysMat2);
 Mat(End_node+1,Inter_node)=1*b/h^2;
 Mat(Inter_node,Inter_node-1:Inter_node+1)=[1 -3 1]*b/h^2;
 Mat(Inter_node,End_node+1)=b/h^2;
 Mat(Inter_node,:)
 size(Mat)
% 
%reaction term
f = @(nn,mm,hh,v) (-1).*(gk.*nn.^4.*(v-ek)+gna.*mm.^3.*hh.*...
    (v-ena)+gl.*(v-el));
for i=1:nT
    %u(1,1:end)=-vstart;
    %Using System matrix
    u(1:nX,i+1)=u(1:nX,i)+(k/c)*(Mat*u(1:nX,i)+...
        f(nn(1:nX,i),mm(1:nX,i),hh(1:nX,i),u(1:nX,i)) );
    %Forward Euler for time dependent ODEs on n,m,h
    nn(1:nX,i+1)=nn(1:nX,i)+k*(an(u(1:nX,i)).*(1-nn(1:nX,i))-bn(u(1:nX,i)).*nn(1:nX,i));
    mm(1:nX,i+1)=mm(1:nX,i)+k*(am(u(1:nX,i)).*(1-mm(1:nX,i))-bm(u(1:nX,i)).*mm(1:nX,i));
    hh(1:nX,i+1)=hh(1:nX,i)+k*(ah(u(1:nX,i)).*(1-hh(1:nX,i))-bh(u(1:nX,i)).*hh(1:nX,i));    
    %u(59,i+1)=-vstart;
    
    u(End_node,i+1)=-vstart;
    u(end,i+1)=-vstart*0;
    u(src_loc,i+1)=-vstart;
    %u(end,i+1)=0;
end


 ymax=max(max(u));
 ymin=min(min(u));
 %figure(1)
 fig=figure('units','normalized','outerposition',[0 0 1 1]);
 vidfile = VideoWriter('MOL_simpleAP_MatrixVersion4.mp4','MPEG-4');
 open(vidfile);

 y1=zeros(1,length(x(1:End_node)));
 %y2=ones(1,length(x(End_node+1:end)))*x(Inter_node);
 y2=linspace(0.5,1,length(x(End_node+1:end)));
 for j=1:50:nT+1
    
    plot3(x(1:End_node),y1,u(1:End_node,j),'-o',...
    y2,x(End_node+1:end)-x(Inter_node+1),u(End_node+1:end,j),'-o','Linewidth',3)
    %legend('Cylinder Part','Branch')
    grid on
    zlim([ymin ymax])
    zlabel('Voltage [mV]')
    xlabel('x pos')
    ylabel('y pos')
    set(gca,'FontSize', 18)
    hold off
    drawnow
    thisFrame = getframe(gcf);
    writeVideo(vidfile, thisFrame);
 end
 
% [Xb,Tb]=meshgrid(x(1:End_node),t);
%  X2=Xb';
%  T2=Tb';
%  
%  [Xb,Tb]=meshgrid(x(End_node+1:end),t);
%  Xb=Xb';
%  Tb=Tb';
%  
%  for j=1:50:nT+1
%      clf
%   subplot(2,2,1)
%   plot(x(1:End_node),u(1:End_node,j),'Linewidth',2);
%   ax = gca;
%   ax.FontSize = 24;
%   xlabel('x pos','fontsize',24);
%   ylabel('V [mV]','fontsize',24);
%   title(sprintf('t = %f [ms]',t(j)),'fontsize',24);
%   ylim([ymin,ymax]);
%   xlim([0,x(End_node)]);
%   
%  
%  subplot(2,2,2)
%  hold on
%  surf(X2,T2,u(1:End_node,:));
%  plot3([x(1) x(End_node)],[t(j) t(j)],[40 40],'r','Linewidth',3);
%  xlim([x(1) x(End_node)])
%  ylim([0 t(end)])
%  hold off
%  view(2)
% shading interp;
% colormap jet;
% ax = gca;
% ax.FontSize = 24;
% hcb=colorbar;
% title(hcb,'voltage [mV]','fontsize',24)
% pos = get(hcb,'Position');
% hcb.Label.Position = [pos(1)/2 pos(2)+1]; % to change its position
% hcb.Label.Rotation = 0; % to rotate the text
% xlabel('x pos','fontsize',24);
% ylabel('time [ms]','fontsize',24);
% zlabel('intensity','fontsize',24);
% set(gca,'FontSize',24)
% 
%  subplot(2,2,3)
%  plot(x(End_node+1:end),u(End_node+1:end,j),'Linewidth',2);
%  ax = gca;
%  ax.FontSize = 24;
%  xlabel('x pos','fontsize',24);
%  ylabel('V [mV]','fontsize',24);
%  title(sprintf('t = %f [ms]',t(j)),'fontsize',24);
%  ylim([ymin,ymax]);
%  xlim([x(End_node+1) x(end)]);
%  
% 
%  subplot(2,2,4)
%  hold on
%  surf(Xb,Tb,u(End_node+1:end,:));
%  plot3([x(End_node+1) x(end)],[t(j) t(j)],[40 40],'r','Linewidth',3);
%  xlim([x(End_node+1) x(end)])
%  ylim([0 t(end)])
%  hold off
%  view(2)
% shading interp;
% colormap jet;
% ax = gca;
% ax.FontSize = 24;
% hcb=colorbar;
% title(hcb,'voltage [mV]','fontsize',24)
% pos = get(hcb,'Position');
% hcb.Label.Position = [pos(1)/2 pos(2)+1]; % to change its position
% hcb.Label.Rotation = 0; % to rotate the text
% xlabel('x pos','fontsize',24);
% ylabel('time [ms]','fontsize',24);
% zlabel('intensity','fontsize',24);
% set(gca,'FontSize',24)

%  drawnow
%  thisFrame = getframe(gcf);
%  writeVideo(vidfile, thisFrame);
% end
end

function out=an(vin)
out=(0.01).*(10-vin)./(exp((10-vin)./10)-1);
end

function out=bn(vin)
out=(0.125).*exp(-vin./80);
end

function out=am(vin)
out=(0.1).*(25-vin)./(exp((25-vin)./10)-1);
end

function out=bm(vin)
out=4.*exp(-vin./18);
end

function out=ah(vin)
out=(0.07).*exp(-vin./20);
end

function out=bh(vin)
out=1./(exp((30-vin)./10)+1);
end

