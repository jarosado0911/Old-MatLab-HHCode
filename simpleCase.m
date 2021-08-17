function [outMM,outHH, outNN,outU,tgrid,xgrid]=simpleCase(n,xi,xf,ti,tf,dt)
%n= number of points on line;

xgrid=linspace(xi,xf,n);
dx = xgrid(2)-xgrid(1);
numTpts=ceil((tf-ti)/dt);
tgrid=linspace(ti,tf,numTpts+1);

a=.01;       %radius
R=40;       %resistance
gk=36;      %potassium ion conductance
ek=-12;     %potassium reversal potential
gna=120;    %sodium ion conductance
ena=112; %115;    %sodium reversal potential
ni=0.5; 
mi=0.4;
hi=0.2;
gl=0.3;     %leak conductance
el=10.6;    %leak potential

b=a/(2*R);
c=0.9e-2;

%b=1;
%c=1;

fprintf('Number of spatial grid points = %i\n',n);
fprintf('Spatial dx = %d\n',dx);
fprintf('Number of temporal grid points = %i\n',numTpts);
fprintf('Temporal dt = %d\n',dt);

%set boundary conditions at one end of line;
u=zeros(n,numTpts); %initialize
nn=zeros(n,numTpts);
mm=zeros(n,numTpts);
hh=zeros(n,numTpts);
nn(1,:)=nn(1,:)+ni;
mm(1,:)=mm(1,:)+mi;
hh(1,:)=hh(1,:)+hi;

vstart=-95;
%for periodic boundary
%for i=1:numTpts
%    u(1,i)=-55e-3*sin(20*tgrid(i));
%end
%i.e. hold -55mV for first half of time
u(1,1:floor(numTpts/8))=-vstart;
%u(1,:)=-vstart;
%u(n,:)=55;
for j=2:n
    for i=1:numTpts
        if j==n
            u(j,i+1)=u(j,i)+(dt/c)*(b/(dx)^2*(0-2*u(j,i)+u(j-1,i))-gk*nn(j,i)^4*(u(j,i)-ek)-gna*mm(j,i)^3*hh(j,i)*(u(j,i)-ena)-gl*(u(j,i)-el));%+I);
        else
            u(j,i+1)=u(j,i)+(dt/c)*(b/(dx)^2*(u(j+1,i)-2*u(j,i)+u(j-1,i))-gk*nn(j,i)^4*(u(j,i)-ek)-gna*mm(j,i)^3*hh(j,i)*(u(j,i)-ena)-gl*(u(j,i)-el));%+I);
        end
        nn(j,i+1)=nn(j,i)+dt*(an(u(j,i))*(1-nn(j,i))-bn(u(j,i))*nn(j,i));
        mm(j,i+1)=mm(j,i)+dt*(am(u(j,i))*(1-mm(j,i))-bm(u(j,i))*mm(j,i));
        hh(j,i+1)=hh(j,i)+dt*(ah(u(j,i))*(1-hh(j,i))-bh(u(j,i))*hh(j,i));
    end
end
outU=u;
outNN=nn;
outMM=mm;
outHH=hh;

ymax=max(max(u));
ymin=min(min(u));
figure(1)
%hold on;
for j=1:n
 plot(tgrid,u(j,:));
 ylim([ymin,ymax]);
 xlim([ti,tf]);
 drawnow
end

figure(2)
[X,Y]=meshgrid(0:1:1,xgrid);
for i=1:200:700000
    Z=repmat(u(:,i),1,2);
    contourf(X,Y,Z);
    drawnow
end

figure(3)
hold on;
for j=1:n
xcoord=ones(1,numTpts+1).*xgrid(j);
plot3(xcoord,tgrid,u(j,:));
view(3);
end

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