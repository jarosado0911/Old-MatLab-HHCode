function [u,mm,nn,hh,x,t]=simpleAP(nX,nT,X,T)
% nX = number of spatial steps along a rod
% nT = number of time steps
% X = Length of rod
% T = end time

%Test case: [u,mm,nn,hh,x,t]=simpleAP(100,200000,1,20);
%Test case: [u,mm,nn,hh,x,t]=simpleAP(200,200000,2,20);
%Test case: [u,mm,nn,hh,x,t]=simpleAP(200,200000,2,30);
if nargin == 0
    [nX,nT, X, T]=deal(120,200000,1,25);
end


k=T/nT; % time step size
h=X/nX; % equidistant spatial discretization width
x=linspace(0,X,nX+1);    % space vector
t=linspace(0,T,nT+1);    % time vector

%----------------------------------------%
%Parameters for AP
a=.001;      %radius
R=10;       %resistance
gk=36;      %potassium ion conductance
ek=-12;     %potassium reversal potential
gna=120;    %sodium ion conductance
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

%Initialize Solutions
u=zeros(nX+1,nT+1);
nn=zeros(nX+1,nT+1);
mm=zeros(nX+1,nT+1);
hh=zeros(nX+1,nT+1);
%set boundary conditions at one end of line;
vstart=-55;
u(1,1:ceil((nT+1)/2))=-vstart;

nn(1,1)=ni;
mm(1,1)=mi;
hh(1,1)=hi;

%Method of Lines Computation
%The j-th row is a position location through time
%The i-th column is a time slice
for i=1:nT
    for j=1:nX
        %Forward Euler for time dependent ODEs on n,m,h
        nn(j,i+1)=nn(j,i)+k*(an(u(j,i))*(1-nn(j,i))-bn(u(j,i))*nn(j,i));
        mm(j,i+1)=mm(j,i)+k*(am(u(j,i))*(1-mm(j,i))-bm(u(j,i))*mm(j,i));
        hh(j,i+1)=hh(j,i)+k*(ah(u(j,i))*(1-hh(j,i))-bh(u(j,i))*hh(j,i));
        
        %If statement is for spatial-temporal
        if j==nX
            %Forward Euler with Centered Difference of 2nd Derivative
            u(j,i+1)=u(j,i)+(k/c)*(b/(h)^2*(0-2*u(j,i)+u(j-1,i))...
                -gk*nn(j,i)^4*(u(j,i)-ek)-gna*mm(j,i)^3*hh(j,i)*...
                (u(j,i)-ena)-gl*(u(j,i)-el));
        elseif j==1
            u(j,i+1)=u(j,i)+(k/c)*(b/(h)^2*(u(j+1,i)-2*u(j,i))...
                -gk*nn(j,i)^4*(u(j,i)-ek)-gna*mm(j,i)^3*hh(j,i)...
                *(u(j,i)-ena)-gl*(u(j,i)-el));
        else
            u(j,i+1)=u(j,i)+(k/c)*(b/(h)^2*(u(j+1,i)-2*u(j,i)+u(j-1,i))...
                -gk*nn(j,i)^4*(u(j,i)-ek)-gna*mm(j,i)^3*hh(j,i)...
                *(u(j,i)-ena)-gl*(u(j,i)-el));
        end
    end
end
u(1,1:ceil((nT+1)/2))=-vstart;

ymax=max(max(u));
ymin=min(min(u));
figure(1)
fig=figure('units','normalized','outerposition',[0 0 1 1]);
%vidfile = VideoWriter('MOL_simpleAP_ifelseVersion.mp4','MPEG-4');
%open(vidfile);
for j=1:500:nT+1
 %subplot(1,2,1)
 plot(x,u(:,j),'Linewidth',2);
 ax = gca;
 ax.FontSize = 24;
 xlabel('x pos','fontsize',24);
 ylabel('V [mV]','fontsize',24);
 title(sprintf('t = %f seconds',t(j)),'fontsize',24);
 ylim([ymin,ymax]);
 xlim([0,X]);
 %pause
%  subplot(1,2,2)
%  plot(x,mm(:,j),x,nn(:,j),x,hh(:,j));
%  xlabel('x pos');
%  ylabel('m,n,h');
%  title(sprintf('t = %f seconds',t(j)));
%  ylim([0,1]);
%  xlim([0,X]);
 drawnow
 %f(j+1)=getframe(fig);
 %writeVideo(vidfile, f(j+1));
end

%Contour Plot with coloring
% figure(2)
% reset(gcf)
% [X,T]=meshgrid(x,t);
% X=X';
% T=T';
% surf(X,T,u);
% ax = gca;
% ax.FontSize = 24;
% view(2)
% shading interp;
% colormap jet;
% hcb=colorbar;
% title(hcb,'voltage [mV]','fontsize',24)
% xlabel('x pos','fontsize',24);
% ylabel('time [ms]','fontsize',24);
% zlabel('voltage [mV]','fontsize',24);

% figure(2)
% surf(X,T,nn);
% shading interp;
% colormap jet;
% colorbar
% xlabel('x pos');
% ylabel('time');
% zlabel('intensity');
% 
% figure(3)
% surf(X,T,mm);
% shading interp;
% colormap jet;
% colorbar
% zlim([0 1])
% xlabel('x pos');
% ylabel('time');
% zlabel('intensity');
% 
% figure(4)
% surf(X,T,hh);
% shading interp;
% colormap jet;
% colorbar
% xlabel('x pos');
% ylabel('time');
% zlabel('intensity');
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

