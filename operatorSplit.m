function [U,x,t]=operatorSplit(nX,nT,x0,xf,t0,tf)

%default
if nargin==0
    [nX, nT,x0, xf, t0, tf]=deal(50,20000,0,1,0,20);
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
a=.01;      % radius
R=40;       % resistance
ni=0.5;     % probability n 
mi=0.4;     % probability m
hi=0.2;     % probability h

% compute coefficients
c=0.009; %c=0.9e-2; %0.009 % capacitance
b=a/(2*R*c);    % diffusion coefficient 
%-------------------------------------------------------------------------%

U=zeros(4*nR,nTime);   % initialize solution

% Initialize one end of rod with fixed voltage
% This is the first row entries of U
vstart=55;
U(1,1:floor(nT/2))=vstart;  % Voltage = 55 at one end
%U(nR,floor(nT/2):floor(0.75*nT))=50;      % Voltage = 0 at other end

% Now initialize m,n,h entries
U(nR+1:2*nR,:)=ni;      %initialize n whole row over time
U(2*nR+1:3*nR,:)=mi;    %initialize m whole row over time
U(3*nR+1:4*nR,:)=hi;    %initialize h whole row over time

% Now operator splitting
for j=1:nTime-1
    
    [Utmp,Ntmp,Mtmp,Htmp,~,~]...
       =reactSolve(U(2:nR-1,j),U(nR+2:2*nR-1,j),U(2*nR+2:3*nR-1,j),U(3*nR+2:4*nR-1,j)...
               ,nX-2,5,x(2),x(nR-1),t(j),t(j)+0.5*k,0, 0);
    U(2:nR-1,j+1)= Utmp(:,end);
    U(nR+2:2*nR-1,j+1)= Ntmp(:,end);
    U(2*nR+2:3*nR-1,j+1)= Mtmp(:,end);
    U(3*nR+2:4*nR-1,j+1)= Htmp(:,end);
    
    %do diffusion solve over time interval with five time steps
    [UUtmp,~,~]=diffusionSolve(U(1,j+1),U(nR,j+1),U(1:nR,j),b,...
                       x0,xf,t(j),t(j+1),nX,5,0);
    U(1:nR,j+1)=UUtmp(1:nR,end); %update
    
    [Utmp,Ntmp,Mtmp,Htmp,~,~]...
       =reactSolve(U(2:nR-1,j+1),U(nR+2:2*nR-1,j+1),U(2*nR+2:3*nR-1,j+1),U(3*nR+2:4*nR-1,j+1)...
               ,nX-2,5,x(2),x(nR-1),t(j),t(j)+0.5*k,0, 0);
    U(2:nR-1,j+1)= Utmp(:,end);
    U(nR+2:2*nR-1,j+1)= Ntmp(:,end);
    U(2*nR+2:3*nR-1,j+1)= Mtmp(:,end);
    U(3*nR+2:4*nR-1,j+1)= Htmp(:,end);
end

figure(1)
[X,T]=meshgrid(x,t);
X=X';
T=T';
surf(X,T,U(1:nR,:));
shading interp;
colormap jet;
colorbar
xlabel('x pos');
ylabel('time');
zlabel('intensity');
end
