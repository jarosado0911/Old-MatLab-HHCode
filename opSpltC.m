function [U,x,t]=opSpltC(nX,nT,x0,xf,t0,tf)

%default
if nargin==0
    [nX, nT,x0, xf, t0, tf]=deal(120,40000,0,1,0,60);
end

k=abs(tf-t0)/nT;         % time step size
h=abs(xf-x0)/nX;         % equidistant spatial discretization width
x=linspace(x0,xf,nX+1);    % space vector
t=linspace(t0,tf,nT+1);    % time vector

fprintf('Number of spatial steps = %i\n',nX);
fprintf('Spatial dx = %d\n',h);
fprintf('Number of time steps = %i\n',nT);
fprintf('Temporal dt = %d\n',k);

nR=length(x);            % size of spatial vector
nTime=nT+1;            % size of time vector

%----------------------------------Parameters-----------------------------%
%Parameters for AP
a=.001;      % radius
R=10;       % resistance

gk=36; %36;      %potassium ion conductance
ek=-12;     %potassium reversal potential
gna=153;   %120 %sodium ion conductance
ena=112; %112; %112;    %115; %sodium reversal potential
ni=0.5; %0.5
mi=0.4;  %0.4
hi=0.2;
gl=0.3;%0.3;     %leak conductance
el=0.6; %0.001;    %10.6 %leak potential

% compute coefficients
c=0.09;%0.005; %c=0.9e-2; %0.009 % capacitance
b=a/(2*R*c); %denom has c    % diffusion coefficient 
%-------------------------------------------------------------------------%

%---For reaction solve-------%
f1 = @(nn,mm,hh,v) -(gk.*(nn.^4).*(v-ek)+gna.*(mm.^3).*hh.*(v-ena)+gl.*(v-el)).*(1/c);
f2 = @(nn,v) an(v).*(1-nn)-bn(v).*nn;
f3 = @(mm,v) am(v).*(1-mm)-bm(v).*mm;
f4 = @(hh,v) ah(v).*(1-hh)-bh(v).*hh;
%-----------------------------%
U=zeros(4*nR,nTime);   % initialize solution

% Initialize one end of rod with fixed voltage
% This is the first row entries of U
vstart=55;
%U(1,1:ceil(nTime/2))=vstart;
U(1,1:end)=0;%vstart;
% Now initialize m,n,h entries
U(nR+1,1)=ni;      %initialize n whole row over time
U(2*nR+1,1)=mi;    %initialize m whole row over time
U(3*nR+1,1)=hi;    %initialize h whole row over time

% Now operator splitting
for j=1:nTime-1
    U(ceil(nX/2),j)=vstart;
    % Diffusion Solve
   [Uout,~,~]=diffSolve(U(1:nR,j),U(1,j+1),U(nR,j+1),b,nX,2,x0,...
                        xf,t(j),t(j+1));
   tmpU=Uout(:,end);
   
    %Now Solve Reaction Equation using forward euler
          U(2:nR-1,j+1)=tmpU(2:nR-1)+k*f1(U(nR+2:2*nR-1,j),...
              U(2*nR+2:3*nR-1,j),U(3*nR+2:4*nR-1,j),tmpU(2:nR-1));
    % Solve n,m,h using forward euler      
       U(nR+1:2*nR,j+1)=U(nR+1:2*nR,j)+k*f2(U(nR+1:2*nR,j),tmpU);
     U(2*nR+1:3*nR,j+1)=U(2*nR+1:3*nR,j)+k*f3(U(2*nR+1:3*nR,j),tmpU);
     U(3*nR+1:4*nR,j+1)=U(3*nR+1:4*nR,j)+k*f4(U(3*nR+1:4*nR,j),tmpU);
end

ymax=max(max(U(1:nR,:)));
ymin=min(min(U(1:nR,:)));

[X2,T2]=meshgrid(x,t);
X2=X2';
T2=T2';

fig=figure('units','normalized','outerposition',[0 0 1 1]);
vidfile = VideoWriter('OpsplitCenter2.mp4','MPEG-4');
open(vidfile);
for k=1:100:nTime
    clf
    subplot(1,2,1)
    plot(x,U(1:nR,k),'Linewidth',2);
    xlabel('x pos','fontsize',24);
    ylabel('V [mV]','fontsize',24);
    title(sprintf('t = %f [ms]',t(k)),'fontsize',24);
    ylim([ymin ymax]);
    %pause
    
    subplot(1,2,2)
    hold on
    surf(X2,T2,U(1:nR,:));
    plot3([x(1) x(end)],[t(k) t(k)],[40 40],'r','Linewidth',3);
    hold off
    view(2);
    shading interp;
    colormap jet;
    hcb=colorbar;
    title(hcb,'voltage [mV]','fontsize',24)
    xlabel('x pos','fontsize',24);
    ylabel('time [ms]','fontsize',24);
    zlabel('voltage [mV]','fontsize',24);
    set(gca,'FontSize',24)
    drawnow
    thisFrame = getframe(gcf);
    writeVideo(vidfile, thisFrame);
end

% figure(2)
% [X,T]=meshgrid(x,t);
% X=X';
% T=T';
% surf(X,T,U(1:nR,:));
% view(2);
% shading interp;
% colormap jet;
% hcb=colorbar;
% title(hcb,'voltage [mV]','fontsize',24)
% xlabel('x pos','fontsize',24);
% ylabel('time [mV]','fontsize',24);
% zlabel('voltage [mV]','fontsize',24);

% figure(2)
% for k=1:1000:nTime-1
%  plot(x,U(1:nR,k));
%  ylim([min(min(U(1:nR,:))) max(max(U(1:nR,:)))]);
%  drawnow
% end

% figure(2)
% surf(X,T,U(nR+1:2*nR,:));
% shading interp;
% colormap jet;
% colorbar
% xlabel('x pos');
% ylabel('time');
% zlabel('intensity');
% 
% figure(3)
% surf(X,T,U(2*nR+1:3*nR,:));
% shading interp;
% colormap jet;
% colorbar
% xlabel('x pos');
% ylabel('time');
% zlabel('intensity');
% 
% figure(4)
% surf(X,T,U(3*nR+1:4*nR,:));
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
