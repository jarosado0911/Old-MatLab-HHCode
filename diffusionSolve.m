function [U,x,t]=diffusionSolve(BC1,BC2,Ui,a,x0,xf,t0,tf,nX,nT,pltOn)
% For this implementation we have BC at the ends of 1d rod
% I use the Crank-Nicolson Method as described in Leveque p. 183
%
% Ui is the initial condition at t0
% a is the diffusion coefficient
% t0 is the start-time
% tf is the end-time
% nT is the number of time steps
% nX is the number of spatial steps
% BC1 and BC2 are the boundary conditions

%Default initialize
if (nargin == 0)
    %note: nX=500, nT=8000 -> instability at (x,t)=(0,0)
    %note: nx=100, nT=8000 -> stable at (x,t) = (0,0)
    
    [nT, nX, a, x0, xf, t0, tf,pltOn]=deal(8000,100,0.015, 0,1,0,100,1);
    BC1=zeros(1,nT+1);
    BC1(1:ceil((nT+1)/2))=50;
    BC2=zeros(1,nT+1);
    BC2(ceil((nT+1)/2):ceil((nT+1)/2)+ceil((nT+1)/4))=-30;
    Ui=zeros(nX+1,1); 
    Ui(1)=BC1(1); 
    Ui(nX+1)=BC2(1);
end

k=abs(tf-t0)/nT;           % time step size
h=abs(xf-x0)/nX;           % equidistant spatial discretization width
x=linspace(x0,xf,nX+1);    % space vector
t=linspace(t0,tf,nT+1);    % time vector

nR=nX+1;                   % number of spatial points                  
nTime=nT+1;                % number of time points

U=zeros(nR,nTime);   % solution output

%initialize
U(:,1)=Ui;
U(1,:)=BC1;
U(nR,:)=BC2;

%parameters and matrices for Crank-Nicolson Solve
r=a*k/(2*h^2);             % constant for tridiagonal matrix
myI=eye(nR-2);               % Identity matrix, careful solve interior equations

%Tridiagonal matrices from Leveque p. 183
A=conv2(myI,[-r (1+2*r) -r],'same');    
B=conv2(myI,[r (1-2*r) r],'same');

%Crank-Nicolson Solve p. 183
for j=1:nT
   %Careful, solve the interior equations
   U(2:nR-1,j+1)=A\(B*U(2:nR-1,j)+r*myI(:,1)*(U(1,j)+U(1,j+1))...
                       +r*myI(:,nR-2)*(U(nR,j)+U(nR,j+1)) );
end

if pltOn==1
ymax=max(max(U));
ymin=min(min(U));

figure(1)
for j=1:50:nTime
 plot(x,U(:,j));
 xlabel('x pos');
 ylabel('u');
 title(sprintf('t = %f seconds',t(j)));
 ylim([ymin,ymax]);
 xlim([x0,xf]);
 drawnow
end

if nT>3
figure(2)
[X,T]=meshgrid(x,t);
X=X';
T=T';
surf(X,T,U);
shading interp;
colormap jet;
colorbar
xlabel('x pos');
ylabel('time');
zlabel('intensity');
end
end
end

