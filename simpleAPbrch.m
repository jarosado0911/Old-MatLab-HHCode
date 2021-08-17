function [u,mm,nn,hh,x,t]=simpleAPbrch(nX,nT,X,T)
% nX = number of spatial steps along a rod
% nT = number of time steps
% X = Length of rod
% T = end time

%Test case: [u,mm,nn,hh,x,t]=simpleAP(100,200000,1,20);
%Test case: [u,mm,nn,hh,x,t]=simpleAP(200,200000,2,20);
%Test case: [u,mm,nn,hh,x,t]=simpleAP(200,200000,2,30);
if nargin == 0
    [nX,nT, X, T]=deal(100,200000,1,25);
end


k=T/nT; % time step size
h=X/nX; % equidistant spatial discretization width
x=linspace(0,X,nX+1);    % space vector
t=linspace(0,T,nT+1);    % time vector

x1=[x(end)+h:h:2*X];

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

%Initialize Solutions
u1=zeros(length(x1),2*nT+1);
nn1=zeros(length(x1),2*nT+1);
mm1=zeros(length(x1),2*nT+1);
hh1=zeros(length(x1),2*nT+1);

%set boundary conditions at one end of line;
vstart=-55;
u(1,1:ceil((nT+1)/2))=-vstart;

nn(1,1)=ni;
mm(1,1)=mi;
hh(1,1)=hi;

% nn(1,:)=nn(1,:)+ni;
% mm(1,:)=mm(1,:)+mi;
% hh(1,:)=hh(1,:)+hi;

%Method of Lines Computation
%The j-th row is a position location through time
%The i-th column is a time slice

for j=1:nX
    for i=1:nT
        %If statement is for spatial-temporal
        if j==nX
            %Forward Euler with Centered Difference of 2nd Derivative
            u(j,i+1)=u(j,i)+(k/c)*(b/(h)^2*(0-2*u(j,i)+u(j-1,i))-gk*nn(j,i)^4*(u(j,i)-ek)-gna*mm(j,i)^3*hh(j,i)*(u(j,i)-ena)-gl*(u(j,i)-el));%+I);
        elseif j==1
            u(j,i+1)=u(j,i)+(k/c)*(b/(h)^2*(u(j+1,i)-2*u(j,i))-gk*nn(j,i)^4*(u(j,i)-ek)-gna*mm(j,i)^3*hh(j,i)*(u(j,i)-ena)-gl*(u(j,i)-el));
        else
            u(j,i+1)=u(j,i)+(k/c)*(b/(h)^2*(u(j+1,i)-2*u(j,i)+u(j-1,i))-gk*nn(j,i)^4*(u(j,i)-ek)-gna*mm(j,i)^3*hh(j,i)*(u(j,i)-ena)-gl*(u(j,i)-el));%+I);
        end
        
        %Forward Euler for time dependent ODEs on n,m,h
        nn(j,i+1)=nn(j,i)+k*(an(u(j,i))*(1-nn(j,i))-bn(u(j,i))*nn(j,i));
        mm(j,i+1)=mm(j,i)+k*(am(u(j,i))*(1-mm(j,i))-bm(u(j,i))*mm(j,i));
        hh(j,i+1)=hh(j,i)+k*(ah(u(j,i))*(1-hh(j,i))-bh(u(j,i))*hh(j,i));
    end
end

ind=50;
u1(1,1:length(u(end-ind,:)))=u(end-ind,:);
nn1(1,1:length(u(end-ind,:)))=nn(end-ind,:);
mm1(1,1:length(u(end-ind,:)))=mm(end-ind,:);
hh1(1,1:length(u(end-ind,:)))=hh(end-ind,:);

for j=2:length(u1(:,1))
    for i=1:length(u1(1,:))-1
        %If statement is for spatial-temporal
        if j==length(u1(:,1))
            %Forward Euler with Centered Difference of 2nd Derivative
            u1(j,i+1)=u1(j,i)+(k/c)*(b/(h)^2*(0-2*u1(j,i)+u1(j-1,i))-gk*nn1(j,i)^4*(u1(j,i)-ek)-gna*mm1(j,i)^3*hh1(j,i)*(u1(j,i)-ena)-gl*(u1(j,i)-el));%+I);
        elseif j==1
            u1(j,i+1)=u1(j,i)+(k/c)*(b/(h)^2*(u1(j+1,i)-2*u1(j,i))-gk*nn1(j,i)^4*(u1(j,i)-ek)-gna*mm1(j,i)^3*hh1(j,i)*(u1(j,i)-ena)-gl*(u1(j,i)-el));
        else
            u1(j,i+1)=u1(j,i)+(k/c)*(b/(h)^2*(u1(j+1,i)-2*u1(j,i)+u1(j-1,i))-gk*nn1(j,i)^4*(u1(j,i)-ek)-gna*mm1(j,i)^3*hh1(j,i)*(u1(j,i)-ena)-gl*(u1(j,i)-el));%+I);
        end
        
        %Forward Euler for time dependent ODEs on n,m,h
        nn1(j,i+1)=nn1(j,i)+k*(an(u1(j,i))*(1-nn1(j,i))-bn(u1(j,i))*nn1(j,i));
        mm1(j,i+1)=mm1(j,i)+k*(am(u1(j,i))*(1-mm1(j,i))-bm(u1(j,i))*mm1(j,i));
        hh1(j,i+1)=hh1(j,i)+k*(ah(u1(j,i))*(1-hh1(j,i))-bh(u1(j,i))*hh1(j,i));
    end
end

ymax=max(max(u));
ymin=min(min(u));
figure(1)
pos1 = get(gcf,'Position'); 
set(gcf,'Position', pos1 - [pos1(3)/2,0,0,0])
for j=1:1000:nT+1
 plot(x,u(:,j));
 xlabel('x pos');
 ylabel('u');
 title(sprintf('t = %f seconds',t(j)));
 ylim([ymin,ymax]);
 xlim([0,X]);
 drawnow
%  subplot(1,2,2)
%  plot(x,mm(:,j),x,nn(:,j),x,hh(:,j));
%  xlabel('x pos');
%  ylabel('m,n,h');
%  title(sprintf('t = %f seconds',t(j)));
%  ylim([0,1]);
%  xlim([0,X]);
%  drawnow
end


figure(2)
set(gcf,'Position', get(gcf,'Position') + [0,0,150,0]); % When Figure(2) is not the same size as Figure(1)
pos2 = get(gcf,'Position');
set(gcf,'Position', pos2 + [pos1(3)/2,0,0,0])
ymax=max(max(u1));
ymin=min(min(u1));
for j=1:1000:length(u1(1,:))
 plot(x1,u1(:,j));
 xlabel('x pos');
 ylabel('u');
 %title(sprintf('t = %f seconds',t(j)));
 ylim([ymin,ymax]);
 xlim([x1(1),x1(end)]);
 drawnow
end

%Contour Plot with coloring
% figure(2)
% [X,T]=meshgrid(x,t);
% X=X';
% T=T';
% surf(X,T,u);
% view(2)
% shading interp;
% colormap jet;
% colorbar
% xlabel('x pos');
% ylabel('time');
% zlabel('intensity');

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

